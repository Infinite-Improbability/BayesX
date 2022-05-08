import argparse
from numpy import arange, sin, cos, full, ceil, savetxt, tile

parser = argparse.ArgumentParser(
    description='Generate mask for given files of ellipses generated by CIAO')
parser.add_argument('maskFiles', nargs='+',
                    help='List of paths to files of ellipses to mask')
parser.add_argument(
    '-s', '--size', help='Size of region being processed, in pixels', default=8192, type=int)
parser.add_argument('-o', '--outputFile', type=argparse.FileType('w'),
                    help='Path to output file', default='mask.txt')
args = parser.parse_args()


class Ellipse:
    def __init__(self, x, y, radius1, radius2, angle):
        super().__init__()
        self.x = x
        self.y = y
        self.radius1 = radius1
        self.radius2 = radius2
        self.angle = angle

    def generate_matrix(self, overdraw=0):
        radius = max(self.radius1, self.radius2) + 1
        r = arange(-radius, radius+1, 1)

        x = tile(r, (len(r), 1))
        y = tile(r, (len(r), 1)).T

        return (x*cos(self.angle) + y*sin(self.angle)) ** 2 / self.radius1 ** 2 + \
            (x*sin(self.angle) - y *
             cos(self.angle)) ** 2 / self.radius2 ** 2 <= (1 + overdraw)

        # return (x_range*cos(self.angle) + y_range*sin(self.angle)) ^ 2 / self.radius1 ^ 2 + \
        #     (x_range*sin(self.angle) - y_range *
        #      cos(self.angle)) ^ 2 / self.radius2 ^ 2 <= 1


def coord2index(x, y, zeroIndex=(0.5,0.5)):
    # Want (0.5,0.5) -> (1,1)
    # Depends on location of origin, etc.
    return [int(max(0, ceil(i+j))) for i,j in zip((x,y), zeroIndex)]

print("Loading ellipses.")
ellipses = []
for filePath in args.maskFiles:
    with open(filePath, 'r') as f:
        for line in f:
            line = line.strip()
            # Currently assuming pixel coordinates. In theory this isn't guaranteed.
            x, y, r1, r2, angle = [
                float(el.strip(' )')) for el in line.split('(')[1].split(',')]
            ellipses.append(Ellipse(x, y, r1, r2, angle))

print("Generating mask array.")
mask = full((args.size,args.size), False)
for el in ellipses:
    elM = el.generate_matrix()
    radius = elM.shape[1] // 2
    x0, y0 = coord2index(el.x-radius, el.y-radius)
    x1, y1 = coord2index(el.x+radius, el.y+radius)
    x1 += 1
    y1 += 1
    
    if el.x-radius < 0:
        elM = elM[-int(el.x-radius)-1:,:]
    if el.x-radius < 0:
        elM = elM[:,-int(el.y-radius)-1:]

    mask[x0:x1, y0:y1] = mask[x0:x1, y0:y1] | elM

print('Exporting array to file.')
savetxt(args.outputFile, mask, fmt='%.1d')
print('File saved as {}.'.format(args.outputFile.name))
