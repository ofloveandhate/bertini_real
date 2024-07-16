import bertini_real as br


br.data.gather_and_save()
surf = br.data.read_most_recent()

br.plot.plot(surf)