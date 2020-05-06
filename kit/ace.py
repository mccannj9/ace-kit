
from kit.contig import Contig


class AceFile(object):
    def __init__(self, filename:str):
        self.filename = filename
        self.file = open(filename)
        as_line = self.file.readline().strip().split()
        self.ncontigs, self.nreads = int(as_line[1]), int(as_line[2])

        # read until first CO line
        self.file.readline()
        self.buffer = [self.file.readline().strip()]

    def __next__(self):
        line = self.file.readline().strip()
        nlines = 1
        while not(line.startswith("CO")):
            self.buffer.append(line)
            line = self.file.readline().strip()
            if nlines > 1:
                last_1, last_2 = self.buffer[-1], self.buffer[-2]
                if not(last_1) and not(last_2):
                    # removes the last extra line
                    self.buffer.pop()
                    break
            nlines += 1
        output = self.buffer
        self.buffer = [line]
        return Contig(output)

    def __repr__(self):
        return f"{self.filename}: {self.ncontigs}, {self.nreads}"

    def readline(self):
        return self.file.readline()

    def close(self):
        self.file.close()
