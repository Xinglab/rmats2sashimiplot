##
## Class for representing figures
##
import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import misopy.sashimi_plot.plot_utils.plot_settings as plot_settings
import misopy.sashimi_plot.plot_utils.plotting as plotting

class Sashimi:
    """
    Representation of a figure.
    """
    def __init__(self, label, output_dir, dimensions=None, png=False,
                 output_filename=None, settings_filename=None,
                 event=None, chrom=None, no_posteriors=False):
        """
        Initialize image settings.
        """
        self.output_ext = ".pdf"
        if png:
            self.output_ext = ".png"
        
        # Plot label, will be used in creating the plot
        # output filename
        self.label = label

        # Set output directory
        self.set_output_dir(output_dir)

        # Plot settings
        self.settings_filename = settings_filename

        if self.settings_filename != None:
            self.settings = plot_settings.parse_plot_settings(settings_filename,
                                                              event=event,
                                                              chrom=chrom,
                                                              no_posteriors=no_posteriors)
        else:
            # Load default settings if no settings filename was given
            self.settings = plot_settings.get_default_settings()

        if output_filename != None:
            # If explicit output filename is given to us, use it
            self.output_filename = output_filename
        else:
            # Otherwise, use the label and the output directory
            self.set_output_filename()
        
        if dimensions != None:
            self.dimensions = dimensions
        else:
            fig_height = self.settings["fig_height"]
            fig_width = self.settings["fig_width"]
            print "Reading dimensions from settings..."
            print " - Height: %.2f" %(float(fig_height))
            print " - Width: %.2f" %(float(fig_width))
            self.dimensions = [fig_width, fig_height]


    def set_output_dir(self, output_dir):
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

    def set_output_filename(self):
        plot_basename = "%s%s" %(self.label, self.output_ext)
        self.output_filename = os.path.join(self.output_dir, plot_basename)

    def setup_figure(self):
        print "Setting up plot using dimensions: ", self.dimensions
        plt.figure(figsize=self.dimensions)

        # If asked, use sans serif fonts
        font_size = self.settings["font_size"]
        if self.settings["sans_serif"]:
            print "Using sans serif fonts."
            plotting.make_sans_serif(font_size=font_size)

    def save_plot(self, plot_label=None):
        """
        Save plot to the output directory. Determine
        the file type.
        """
        if self.output_filename == None:
            raise Exception, "sashimi_plot does not know where to save the plot."
        output_fname = None
        if plot_label is not None:
            # Use custom plot label if given
            ext = self.output_filename.rsplit(".")[0]
            dirname = os.path.dirname(self.output_filename)
            output_fname = \
                os.path.dirname(dirname, "%s.%s" %(plot_label, ext))
        else:
            output_fname = self.output_filename
        print "Saving plot to: %s" %(output_fname)
        plt.savefig(output_fname)
            
        
