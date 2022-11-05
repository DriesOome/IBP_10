def read_fasta(fp):
  print(">Reading fasta file")
  fastaDict = {}
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: fastaDict[name] = ''.join(seq)
      name, seq = line, []
    else:
      seq.append(line)
  if name: fastaDict[name] = ''.join(seq)
  print(">DONE")
  return fastaDict
