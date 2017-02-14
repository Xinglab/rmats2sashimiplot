#!/usr/bin/python
import os
import sashimi_plot
from subprocess import call

sub_path = 'sashimi_plot_new/sashimi_plot/'
path = os.path.abspath('sashimi_plot_new/sashimi_plot')
call(['cd',path], shell = True)
call(['index_gff','--index' ,sub_path + 'test-data/events.gff', sub_path + 'test-data/event-data/'])

call(['sashimi_plot', '--plot','-event' ,'\"chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-\"' ,'test-data/event-data/', 'settings/sashimi_plot_settings.txt', '--output','-dir', 'test-plot'], shell = True)
