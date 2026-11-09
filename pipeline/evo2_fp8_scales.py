"""Use the FP8 scales stored in an Evo 2 checkpoint (1B / 20B / 40B).

Evo 2's Hyena blocks compute their input projection as an FP8 GEMM through Transformer Engine. TE keeps the per-tensor
scales it was trained with in `blocks.N.projections._extra_state`, but vortex's loader reports those entries as unused
("Extra keys in state_dict") and TE then derives scales of its own. `patch(model, checkpoint)` replaces each such
projection with a module that casts input and weight to E4M3 with the stored scales and multiplies them with
torch._scaled_mm (native FP8 on Hopper, e.g. H100).

TESTED, NOT NEEDED ON HOPPER (2026-11-09, 40B on 2 x H100, Fir): the stock TE path already gives evo2's reference
(loss 0.2160 / 91.67%); with the stored scales and an fp32 GEMM the loss is 1.4042 (mmaemu reports 1.4023 for the 40B
on non-Hopper FP8), and the residual-stream jump to ~1e16 at block 23 is there either way, so it belongs to the model.
Use this only on GPUs whose FP8 GEMMs differ from Hopper's, and then with mmaemu's emulated GEMM, not the fp32 one.

Adapted from mmaemu's examples/evo2_offload.py (Xun Ding and Ziyuan Qin, Apache-2.0,
https://github.com/astding/mmaemu): read_te_scales and the Projection module's "native" mode.
"""
import io
import torch
import torch.nn.functional as F

E4M3_MAX = 448.0


def read_te_scales(checkpoint):
    """{block index: FP8 scale [input, weight]} from the TE `projections._extra_state` entries of the checkpoint."""
    torch.serialization.add_safe_globals([io.BytesIO])
    state = torch.load(checkpoint, map_location="cpu", mmap=True, weights_only=True)
    scales = {}
    for key, extra in state.items():
        if not key.endswith(".projections._extra_state") or extra is None:
            continue
        if isinstance(extra, io.BytesIO):
            extra.seek(0); meta = torch.load(extra, map_location="cpu", weights_only=True)
        elif isinstance(extra, torch.Tensor):    # newer TE versions serialise the metadata into a uint8 tensor
            meta = torch.load(io.BytesIO(extra.numpy().tobytes()), map_location="cpu", weights_only=True)
        else:
            meta = extra
        scales[int(key.split(".")[1])] = meta["scale_fwd"][:2].float().clone()
    return scales


class StoredScaleFP8(torch.nn.Module):
    """A Hyena block's input projection as TE computes it, with the checkpoint's stored scales."""
    dequant = False   # set once torch._scaled_mm has been refused; then every projection takes the fp32 path

    def __init__(self, weight, scale):
        super().__init__()
        self.weight = weight                     # [N, K] bf16 parameter, on its block's device
        self.register_buffer("scale", scale.to(weight.device), persistent=False)

    def forward(self, x):
        s = self.scale.to(x.device); s_inv = torch.reciprocal(s)
        wq = (self.weight.detach().float() * s[1]).clamp(-E4M3_MAX, E4M3_MAX).to(torch.float8_e4m3fn)
        shape, L = x.shape, x.shape[-2]
        pad = (16 - L % 16) % 16                 # vortex's TE path pads the sequence length to a multiple of 16
        if pad:
            x = F.pad(x, (0, 0, 0, pad))
        xq = (x.float() * s[0]).clamp(-E4M3_MAX, E4M3_MAX).to(torch.float8_e4m3fn).reshape(-1, x.shape[-1])
        try:
            if StoredScaleFP8.dequant: raise RuntimeError("dequant mode")
            y = torch._scaled_mm(xq, wq.t(), scale_a=s_inv[0], scale_b=s_inv[1], out_dtype=torch.bfloat16,
                                 use_fast_accum=True)
        except RuntimeError:
            # cuBLASLt in some builds rejects these FP8 shapes (CUBLAS_STATUS_NOT_SUPPORTED): multiply the same E4M3
            # values in fp32 instead (same quantisation and scales, fp32 rather than fast accumulation)
            StoredScaleFP8.dequant = True
            y = ((xq.float() @ wq.float().t()) * (s_inv[0] * s_inv[1])).to(torch.bfloat16)
        return y.reshape(*shape[:-2], -1, self.weight.shape[0])[..., :L, :], None   # TELinear's (output, bias)


def patch(model, checkpoint):
    """Replace every FP8 projection of a loaded vortex StripedHyena (evo2's model.model) in place; returns the count."""
    scales = read_te_scales(checkpoint)
    for i, s in scales.items():
        p = model.blocks[i].projections
        b = getattr(p, "bias", None)               # TE keeps an empty tensor for bias=False
        assert b is None or b.numel() == 0 or not b.abs().any(), f"block {i}: Evo 2's FP8 projections have no bias"
        model.blocks[i].projections = StoredScaleFP8(p.weight, s)
    return len(scales)
