#--- Documentation block -------------------------------------------------------
#
#  convert WMO standard BUFR to NCEP BUFR for only GPS RO data  
#
#  Note: 
#    ##  bufr.write_subset(hdr[38:39],'GPHTST',end=True)
#    this python script do NOT pass out 'GPHTST' to the NCEP BUFR, althrough WMO BUFR 
#    has value but the missing values. Because of the bug difficulty taken in count.
#    Wait for future debugging.
#
#  Created by Deng-Shun Chen
#  
#  Log : 
#    2018-05-10  Deng-Shun  Created
# 
#
# CopyRight @ Deng-Shun Chen
#
#--- End of Documentation block ------------------------------------------------
#!/usr/bin/env python
from __future__ import print_function
import ncepbufr
from pybufrkit.decoder import Decoder
from pybufrkit.renderer import FlatTextRenderer, NestedTextRenderer, FlatJsonRenderer, NestedJsonRenderer
import numpy as np
import json 
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

#---------------------------------------------------------------------------------#
parser = ArgumentParser(description = 'Convert WMO BUFR to NCEP BUFR',formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile',help='Input WMO BUFR filename',type=str,required=True)
parser.add_argument('-o','--outfile',help='Output NCEP BUFR filename',type=str,required=False,default=None)
parser.add_argument('-t','--bufrtable',help='NCEP BUFR table file',type=str,required=False,default='gpsbufr.table')
parser.add_argument('-a','--append',help='append to existing NCEP BUFR file',action='store_true',required=False)

args = parser.parse_args()

infilename  = args.infile
outfilename = args.infile+'_ncepbufr' if args.outfile is None else args.outfile
bufrtable   = args.bufrtable
append      = args.append

print('Converting ',infilename,' to ',outfilename)

_nmaxseq = 255 # max size of sequence in message

#---------------------------------------------------------------------------------#
# Open files
#---------------------------------------------------------------------------------#
#** Open WMO bufr file 
decoder = Decoder()
with open(infilename,'rb') as ins:
  bufr_message = decoder.process(s=ins.read(),file_path=infilename)

#** open NCEP bufr file
if append:
  bufr = ncepbufr.open(outfilename,'a')
else:
  bufr = ncepbufr.open(outfilename,'w',table=bufrtable)

#---------------------------------------------------------------------------------#
# Read WMO Standard BUFR
#---------------------------------------------------------------------------------#
#** get whole data
data = FlatJsonRenderer().render(bufr_message)[3][2][0]

#** get header
hdr = bufr.missing_value*np.ones((44),np.float)
hdr[0:37] = np.asarray(data[0:37])
yyyymmddhh ='%04i%02i%02i%02i%02i' % tuple(hdr[6:11])

#** get number of replications of ROSEQ1
nreps_this_ROSEQ1 = np.asarray(data[37:38])

#** get data1b (roseq1)
data1b = bufr.missing_value*np.ones(([_nmaxseq,int(nreps_this_ROSEQ1)]),np.float)
nreps_this_ROSEQ2 = bufr.missing_value*np.ones(([int(nreps_this_ROSEQ1)]),np.float)

n = 38  # skip header (should be no change)
for k in range(int(nreps_this_ROSEQ1)):
  for x in range(3):
    ii = n + x
    jj = ii + 1
    data1b[x,k] = np.asarray(data[ii:jj])
  ii = ii + 1
  jj = jj + 1 
  nreps_this_ROSEQ2[k] = np.asarray(data[ii:jj])

  l = n+4
  for i in range(int(nreps_this_ROSEQ2[k])):
    m = 6*(i+1)-3
    j = 6*i + l
    for x in range(6): 
      ii = j+x 
      jj = ii + 1
      mx = m+x
      if None in data[ii:jj]:
        data1b[mx,k] = bufr.missing_value
      else:
        data1b[mx,k] = np.asarray(data[ii:jj])

  ii = ii + 1
  jj = jj + 1 
  data1b[mx+1,k] = np.asarray(data[(ii+1):(jj+1)])
  n = n + 5 + int(nreps_this_ROSEQ2[k])*6 

## get data2a (roseq3)
nreps_this_ROSEQ3 = np.asarray(data[(n):(n+1)])
data2a = bufr.missing_value*np.ones(([_nmaxseq,int(nreps_this_ROSEQ3)]),np.float)

n = n + 1
for k in range(int(nreps_this_ROSEQ3)):
  for x in range(6):
    ii = n + x
    jj = ii + 1
    if None in data[ii:jj]:
      data2a[x,k] = bufr.missing_value
    else:
      data2a[x,k] = np.asarray(data[ii:jj])
  n = n + 6

## get data2b (roseq4)
nreps_this_ROSEQ4 = np.asarray(data[(n):(n+1)])
data2b = bufr.missing_value*np.ones(([_nmaxseq,int(nreps_this_ROSEQ4)]),np.float)

n = n + 1
for k in range(int(nreps_this_ROSEQ4)):
  for x in range(10):
    ii = n + x
    jj = ii + 1
#    print('ii=',ii,'jj=',jj,'data=',data[ii:jj])
    if None in data[ii:jj]:
      data2b[x,k] = bufr.missing_value
    else:
      data2b[x,k] = np.asarray(data[ii:jj])
  n = n + 10

hdr[37:44] = np.asarray([bufr.missing_value if v is None else v for v in data[n:n+7]])

#---------------------------------------------------------------------------------#
# Write out NCEP BUFR
#---------------------------------------------------------------------------------#
#** open message
idate=yyyymmddhh[0:10] # cycle time: YYYYMMDDHH
subset='NC003010'

print('subset/idate = ',subset,idate)
print('nreps_this_ROSEQ1 = ',int(nreps_this_ROSEQ1))
print('nreps_this_ROSEQ3 = ',int(nreps_this_ROSEQ3))
print('nreps_this_ROSEQ4 = ',int(nreps_this_ROSEQ4))

# open message
bufr.open_message(subset, idate)

#** write subset **
# write seq 0
bufr.write_subset(hdr[0:4],'SIIDSEQ',seq=True)
bufr.write_subset(hdr[4:6],'SWID TSIG')
bufr.write_subset(hdr[6:9],'YYMMDD',seq=True)
bufr.write_subset(hdr[9:11],'HHMM',seq=True)
bufr.write_subset(hdr[11:13],'SECO QFRO')
bufr.write_subset(hdr[20:22],'SCLF PTID')
bufr.write_subset(hdr[28:29],'TISE')
bufr.write_subset(hdr[29:31],'LTLONH',seq=True)
bufr.write_subset(hdr[34:37],'ELRC BEARAZ GEODU')
bufr.write_subset(hdr[37:38],'VSAT')

locplat = bufr.missing_value*np.ones(([3,3]),np.float)
locplat[:,0] = hdr[14:17]
locplat[:,1] = hdr[22:25]
locplat[:,2] = hdr[31:34] 
bufr.write_subset(locplat,'LOCPLAT',seq=True)

spdplat = bufr.missing_value*np.ones(([3,2]),np.float)
spdplat[:,0] = hdr[17:20]
spdplat[:,1] = hdr[25:28]
bufr.write_subset(spdplat,'SPDPLAT',seq=True)

tmp = bufr.missing_value*np.ones(([2,2]),np.float)
tmp[:,0] = hdr[39:41]
tmp[:,1] = hdr[41:43]
bufr.write_subset(tmp,'PRES FOST',rep=True)

tmp = bufr.missing_value*np.ones(([1,2]),np.float)
tmp[0,0] = hdr[13:14]
tmp[0,1] = hdr[43:44]
bufr.write_subset(tmp,'PCCF',rep=True)

#bufr.write_subset(hdr[38:39],'GPHTST')

#**  make a template using drfini
bufr.drfini([nreps_this_ROSEQ1],'(ROSEQ1)')
bufr.drfini(nreps_this_ROSEQ2,'{ROSEQ2}')
bufr.drfini([nreps_this_ROSEQ3],'(ROSEQ3)')
bufr.drfini([nreps_this_ROSEQ4],'(ROSEQ4)')

# write drf ROSEQ1 (L1B)
bufr.write_subset(data1b,'ROSEQ1',seq=True)

# write drf ROSEQ3 (L2A)
bufr.write_subset(data2a,'ROSEQ3',seq=True)

# write drf ROSEQ4 (L2B)
bufr.write_subset(data2b,'ROSEQ4',seq=True,end=True)

#** close message
bufr.close_message()

#** close ncepbufr file
bufr.close()

print('Successful End')

