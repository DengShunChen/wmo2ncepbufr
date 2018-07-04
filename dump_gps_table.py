from __future__ import print_function
import ncepbufr

# read gpsro file.

bufr = ncepbufr.open('gpsbufr')
bufr.dump_table('gpsbufr.table')

