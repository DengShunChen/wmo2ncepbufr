import ncepbufr

bufr = ncepbufr.open('amsuabufr')
bufr.dump_table('amsuabufr.table')
bufr.close()

