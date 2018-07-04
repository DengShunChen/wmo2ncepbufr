from __future__ import print_function
import ncepbufr
import numpy as np

hdrstr ='YEAR MNTH DAYS HOUR MINU PCCF ELRC SAID PTID GEODU'
roseq2str ='CLATH CLONH BEARAZ PCCF'

# read gpsro file.

bufr = ncepbufr.open('gpsbufr')
bufr.print_table()
while bufr.advance() == 0:
    print(bufr.msg_counter, bufr.msg_type, bufr.msg_date)
    while bufr.load_subset() == 0:
        hdr = bufr.read_subset(hdrstr).squeeze()
        yyyymmddhh ='%04i%02i%02i%02i%02i' % tuple(hdr[0:5])
        satid = int(hdr[7])
#        if satid == 740:
        ptid = int(hdr[8])
        nreps_this_ROSEQ2 = bufr.read_subset('{ROSEQ2}').squeeze()
        nreps_this_ROSEQ1 = len(nreps_this_ROSEQ2)
        print(nreps_this_ROSEQ1,nreps_this_ROSEQ2)
        data1b = bufr.read_subset('ROSEQ1',seq=True) # bending angle
        data2a = bufr.read_subset('ROSEQ3',seq=True) # refractivity
        levs_bend = data1b.shape[1]
        levs_ref = data2a.shape[1]
        if levs_ref != levs_bend:
            print('skip report due to bending angle/refractivity mismatch')
            continue
        print('sat id,platform transitter id, levels, yyyymmddhhmm =',\
        satid,ptid,levs_ref,yyyymmddhh)
        print('k, lat, lon,  height, ref, bend:')
        for k in range(int(nreps_this_ROSEQ1)):
          rlat = data1b[0,k]
          rlon = data1b[1,k]
          height = data2a[0,k]
          ref = data2a[1,k]
          for i in range(int(nreps_this_ROSEQ2[k])):
              m = 6*(i+1)-3
              freq = data1b[m,k]
              bend = data1b[m+2,k]
              # look for zero frequency bending angle ob
              if int(freq) == 0: break
          print(k,rlat,rlon,height,ref,bend)
    # only loop over first 6 subsets
    if bufr.msg_counter == 5: break
bufr.close()

#print(np.asarray(data1b[0:2,:])) 
print(data1b.shape) 

bufr2 = ncepbufr.open('gpsbufr_new','w',table='gpsbufr.table')
idate=yyyymmddhh[0:10] # cycle time: YYYYMMDDHH
subset='NC003010'
bufr2.open_message(subset, idate)

# initialize drf 
print(int(nreps_this_ROSEQ2[0]))
bufr2.drfini([nreps_this_ROSEQ1],'(ROSEQ1)')
bufr2.drfini([nreps_this_ROSEQ1],'(ROSEQ3)')
bufr2.drfini(nreps_this_ROSEQ2,'{ROSEQ2}')

# write seq 0
bufr2.write_subset(hdr,hdrstr)
bufr2.write_subset(data1b,'ROSEQ1',seq=True)
bufr2.write_subset(data2a,'ROSEQ3',seq=True,end=True)

bufr2.close_message()
bufr2.close()

