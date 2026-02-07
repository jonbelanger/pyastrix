import numpy as np
from astropy.visualization import ImageNormalize, AsinhStretch

class ImageNormalizer:
    def __init__(self, mode="auto", gamma=1.1):
        """
        mode: 'auto', 'linear', 'stretched'
        """
        self.mode = mode
        self.gamma = gamma

    def _looks_linear(self, data):
        flat = data.flatten().astype(np.float32)
        p1, p50, p99 = np.percentile(flat, [1, 50, 99])
        # linear data has huge dynamic range and low midtone occupancy
        return (p99 / max(p1, 1e-6)) > 1e3 and p50 < 0.05 * p99

    def normalize(self, data):
        if data.ndim == 3:
            rgb = np.moveaxis(data, 0, -1).astype(np.float32)
        else:
            rgb = np.stack([data]*3, axis=-1).astype(np.float32)

        mode = self.mode
        if mode == "auto":
            mode = "linear" if self._looks_linear(rgb) else "stretched"

        if mode == "linear":
            return self._normalize_linear(rgb)
        else:
            return self._normalize_stretched(rgb)

    def _normalize_linear(self, rgb):
        out = np.zeros_like(rgb, dtype=np.float32)
        stretch = AsinhStretch(a=0.01)

        for c in range(3):
            ch = rgb[:, :, c]
            vmin = np.percentile(ch, 0.25)
            vmax = np.percentile(ch, 99.7)

            norm = ImageNormalize(
                ch,
                vmin=vmin,
                vmax=vmax,
                stretch=stretch
            )
            out[:, :, c] = norm(ch)

        return np.clip(out, 0, 1)

    def _normalize_stretched(self, rgb):
        lo = np.percentile(rgb, 0.5)
        hi = np.percentile(rgb, 99.5)
        rgb = (rgb - lo) / (hi - lo)
        rgb = np.clip(rgb, 0, 1)
        rgb = np.power(rgb, 1 / self.gamma)
        return rgb
