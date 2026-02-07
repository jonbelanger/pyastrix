import numpy as np
from adjustText import adjust_text
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class ImageWriter:
    def __init__(self, dpi=200):
        self.dpi = dpi
        self.ny = None
        self.nx = None

    def begin(self, rgb):
        """
        rgb: normalized image array (ny, nx, 3)
        """
        self.ny, self.nx = rgb.shape[:2]

        dpi = self.dpi
        fig = plt.figure(figsize=(self.nx / dpi, self.ny / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])

        # force the image to exactly original pixel size
        ax.imshow(rgb, origin="lower", extent=(0, self.nx, 0, self.ny))
        ax.set_xlim(0, self.nx)
        ax.set_ylim(0, self.ny)
        ax.set_aspect('equal')
        ax.axis('off')

        return fig, ax

    def finish(self, fig, ax, filename, texts=None):
        # run adjustText if there are labels
        if texts:
            adjust_text(
                texts,
                only_move={'points':'y','text':'y'},
                arrowprops=dict(arrowstyle='-', color='green', lw=0.5),
                ax=ax
            )

        # draw the canvas
        fig.canvas.draw()

        # grab RGBA buffer
        w, h = fig.canvas.get_width_height()
        buf = np.frombuffer(fig.canvas.renderer.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
        rgb = buf[:, :, :3]

        # crop to original image dimensions
        cropped = rgb[:self.ny, :self.nx, :]

        # save
        plt.imsave(filename, cropped)
        plt.close(fig)