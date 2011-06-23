def process_commands(text):
    "Scan text for animation commands and process any found."
    chunks=text.split("ANIMATION(")
    results=chunks[0]
    for chunk in chunks[1:]:
        pieces=chunk.split(")")
        f=pieces[0]+".Rnw"
        result=process_template(f)
        results += result+")".join(pieces[1:])
    return results
def process_template(f):
    "Scan template file for %s and fill in."
    lines=open(f).readlines()
    items=lines[0].replace("%","").split()
    text="".join(lines[1:])
    slides=[text.replace("%s",i) for i in items]
    return '\n'.join(slides)
if __name__ == "__main__":
    commands = open("2011-07-01-HOCKING-icml-clusterpath-slides.Ranim").read()
    f=open("2011-07-01-HOCKING-icml-clusterpath-slides.Rnw","w")
    f.write(process_commands(commands))
