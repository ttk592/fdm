e.g. to generate a video use mencoder:

 $ mencoder "mf://out/*.png" -mf fps=25 -o out.avi -ovc lavc -lavcopts vcodec=mpeg4
