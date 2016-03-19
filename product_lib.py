import matplotlib
matplotlib.use('Agg')
import pylab as pl


def produceQuickImagePlot(img, title, x_label, y_label):
    #return;
    pl.figure('quick image plot', facecolor='white', figsize=(8, 6))
    pl.cla()
    pl.title(title, fontsize=14)
    pl.xlabel(x_label, fontsize=12)
    pl.ylabel(y_label, fontsize=12)
    #pl.imshow(img, aspect='auto')
    pl.imshow(img, vmin=0, vmax=256, aspect='auto', cmap="gray")
#     pl.colorbar()
    #pl.set_cmap('spectral')
    pl.show()
