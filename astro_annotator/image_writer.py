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
                only_move={'points':'y','text':'xy'},
                expand_text=(1.2, 1.2),
                expand_points=(1.2, 1.2),
                arrowprops=dict(arrowstyle='-', lw=0, color='none'),  # suppress adjust_text arrows
                ax=ax
            )
            
            # Draw manual connector lines from text to just outside circle
            for txt in texts:
                circle_center = getattr(txt, 'circle_center', None)
                if circle_center is None:
                    continue
                
                # Get text position in data coordinates
                text_x, text_y = txt.get_position()
                cx, cy = circle_center
                
                # Calculate direction from circle to text
                dx = text_x - cx
                dy = text_y - cy
                dist = np.sqrt(dx**2 + dy**2)
                
                if dist > 0:
                    # Normalize and move endpoint a few pixels from circle edge
                    # Circle marker has radius ~2.5 (markersize=5), add ~3-4 pixels
                    offset = 5.5  # pixels from circle center to line endpoint
                    unit_x = dx / dist
                    unit_y = dy / dist
                    line_end_x = cx + unit_x * offset
                    line_end_y = cy + unit_y * offset
                else:
                    line_end_x = cx
                    line_end_y = cy
                
                # Draw line from text to just outside circle
                ax.plot(
                    [text_x, line_end_x],
                    [text_y, line_end_y],
                    color='green',
                    linewidth=0.7,
                    linestyle='-',
                    alpha=0.7,
                    zorder=0  # draw behind text
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