template = open("2011-07-01-HOCKING-icml-clusterpath-slides-in.Rnw").read()
slide_tmp = open("panels.Rnw").read()
slides = [slide_tmp%(s,s) for s in [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]]
slides_sweave = '\n'.join(slides)
sweave = template.replace("COMPARE_NORM_WEIGHTS",slides_sweave)
f=open("2011-07-01-HOCKING-icml-clusterpath-slides.Rnw","w")
f.write(sweave)
