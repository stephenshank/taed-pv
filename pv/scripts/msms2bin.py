#!/usr/bin/env python
"""
Converts the MSMS .vert/.face pair into a compact binary format optimized for
file size.
"""
import sys
import os
from struct import pack 

def print_usage(code=0):
    print 'usage: msm2bin.py <input> [output]'
    print '  input may either point to the .vert or .face file generated by msms'
    print '  when the output filename is omitted, it defaults to input.bin with'
    print '  the file extension removed, e.g. surface.vert -> surface.bin'
    sys.exit(code);

def exchange_extension(filename, new_ext):
    return '%s.%s' % (os.path.splitext(filename)[0], new_ext)


def check_exists(filename):
    if not os.path.exists(filename):
        print '"%s" does not exists' % filename
        return False
    return True

def check_readable(filename):
    if not check_exists(filename):
        return False
    if not os.access(filename, os.R_OK):
        print '"%s" is not readable' % filename
        return False
    return True

def read_and_pack_vert(filename, out_fd):
    with open(filename, 'rb') as verts:
        header = True
        data = bytes()
        for line in verts:
            line = line.strip()
            if line.startswith('#'):
                continue
            if len(line) == 0:
                continue
            if header:
                num_verts = int(line.split()[0])
                out_fd.write(pack('!I', num_verts))
                header = False
                continue
            x, y, z, nx, ny, nz = line.split()[:6]
            data += pack('!ffffff', float(x), float(y), float(z), 
                         float(nx), float(ny), float(nz))
        out_fd.write(data)
def read_and_pack_faces(filename, out_fd):
    with open(filename, 'rb') as faces:
        header = True
        data = bytes()
        for line in faces:
            line = line.strip()
            if line.startswith('#'):
                continue
            if len(line) == 0:
                continue
            if header:
                num_faces = int(line.split()[0])
                out_fd.write(pack('!I', num_faces))
                header = False
                continue
            idx0, idx1, idx2  = line.split()[:3]
            data += pack('!III', int(idx0), int(idx1), int(idx2))
        out_fd.write(data)



def main(args):
    BIN_VERSION = 1
    if len(args) < 2:
        print_usage(-1)
    vert_filename = exchange_extension(args[1], 'vert')
    face_filename = exchange_extension(args[1], 'face')
    output_filename = exchange_extension(args[1], 'bin')
    if len(args) == 3:
        output_filename = args[2]
    if not check_readable(vert_filename) or not check_readable(face_filename):
        sys.exit(-1)
    with open(output_filename, 'wb') as out_file:
        out_file = open(output_filename, 'wb')
        out_file.write(pack('!I', BIN_VERSION))
        read_and_pack_vert(vert_filename, out_file)
        read_and_pack_faces(face_filename, out_file)

if __name__ == '__main__':
    main(sys.argv)


