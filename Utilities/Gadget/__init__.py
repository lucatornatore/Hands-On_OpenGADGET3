from __future__ import print_function
# routines for reading headers and data blocks from Gadget snapshot files
# usage e.g.:
#
# import Gadget as G
# header = G.snapshot_header("snap_063.0") # reads snapshot header
# print header.massarr
# mass = rs.read_block("snap_063","MASS",parttype=5) # reads mass for particles of type 5, using block names should work for both format 1 and 2 snapshots
# print "mass for", mass.size, "particles read"
# print mass[0:10]
#
# before using read_block, make sure that the description (and order if using format 1 snapshot files) of the data blocks
# is correct for your configuration of Gadget 
#
# for mutliple file snapshots give e.g. the filename "snap_063" rather than "snap_063.0" to read_block
# for snapshot_header the file number should be included, e.g."snap_063.0", as the headers of the files differ
#
# the returned data block is ordered by particle species even when read from a multiple file snapshot

import numpy as np
import os
import sys
import math
  
# ----- class for snapshot header ----- 



class snapshot_header:
  def __init__(self, filename):

    if os.path.exists(filename):
      curfilename = filename
    elif os.path.exists(filename+".0"):
      curfilename = filename+".0"
    else:
      raise ValueError("File not found:"+filename)
      #sys.exit()
      
    self.filename = filename  
    f = open(curfilename,'rb')    
    blocksize = np.fromfile(f,dtype=np.int32,count=1)
    if blocksize[0] == 8:
      swap = 0
      format = 2
    elif blocksize[0] == 256:
      swap = 0
      format = 1  
    else:
      blocksize.byteswap(True)
      if blocksize[0] == 8:
        swap = 1
        format = 2
      elif blocksize[0] == 256:
        swap = 1
        format = 1
      else:
        raise ValueError("incorrect file format encountered when reading header of"+ filename)
    
    self.format = format
    self.swap = swap
    
    if format==2:
      f.seek(16, os.SEEK_CUR)
    
    self.npart = np.fromfile(f,dtype=np.int32,count=6)
    self.massarr = np.fromfile(f,dtype=np.float64,count=6)
    self.time = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.redshift = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.sfr = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.feedback = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.nall = np.fromfile(f,dtype=np.int32,count=6)
    self.cooling = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.filenum = (np.fromfile(f,dtype=np.int32,count=1))[0]
    self.boxsize = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_m = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.omega_l = (np.fromfile(f,dtype=np.float64,count=1))[0]
    self.hubble = (np.fromfile(f,dtype=np.float64,count=1))[0]
    
    if swap:
      self.npart.byteswap(True)
      self.massarr.byteswap(True)
      self.time = self.time.byteswap()
      self.redshift = self.redshift.byteswap()
      self.sfr = self.sfr.byteswap()
      self.feedback = self.feedback.byteswap()
      self.nall.byteswap(True)
      self.cooling = self.cooling.byteswap()
      self.filenum = self.filenum.byteswap()
      self.boxsize = self.boxsize.byteswap()
      self.omega_m = self.omega_m.byteswap()
      self.omega_l = self.omega_l.byteswap()
      self.hubble = self.hubble.byteswap()
     
    f.close()
 
# ----- find offset and size of data block ----- 

def find_block(filename, format, swap, block, block_num, only_list_blocks=False):
  if (not os.path.exists(filename)):
      raise ValueError("File not found:"+filename)

  f = open(filename,'rb')
  f.seek(0, os.SEEK_END)
  filesize = f.tell()
  f.seek(0, os.SEEK_SET)
  
  found = False
  curblock_num = 1
  while ((not found) and (f.tell()<filesize)):
    if format==2:
      f.seek(4, os.SEEK_CUR)
      curblock = f.read(4)
      curblock=curblock.decode('utf-8')
      if (block == curblock):
        found = True
      f.seek(8, os.SEEK_CUR)  
    else:
      if curblock_num==block_num:
        found = True


    curblocksize = (np.fromfile(f,dtype=np.uint32,count=1))[0]
    if swap:
      curblocksize = curblocksize.byteswap()
    
    # - print some debug info about found data blocks -
    #if format==2:
    #  print curblock, curblock_num, curblocksize
    #else:
    #  print curblock_num, curblocksize
    
    if only_list_blocks:
      if format==2:
        print(curblock_num,curblock,f.tell(),curblocksize)
      else:
        print(curblock_num,f.tell(),curblocksize)
      found = False
        
    
    if found:
      blocksize = curblocksize
      offset = f.tell()
    else:
      f.seek(curblocksize, os.SEEK_CUR)
      blocksize_check = (np.fromfile(f,dtype=np.uint32,count=1))[0]
      if swap: blocksize_check = blocksize_check.byteswap()
      if (curblocksize != blocksize_check):
        raise ValueError("something wrong")
      curblock_num += 1
  f.close()
      
  if ((not found) and (not only_list_blocks)):
    raise ValueError("block not found")

  if (not only_list_blocks):
    return offset,blocksize
 
# ----- read data block -----
 
def read_block(filename, blockd, parttype=-1, physical_velocities=False, arepo=0, no_masses=False, verbose=False, muppi=1, nMetals = 15):
  """
       USE: read_block('snapfile','BLCK')

       AUXILIARY:
       parttype=-1                only load GadgetType n (-1 = all)
       physical_velocities=True   multiply velocities by sqrt(expansion_factor)
       nMetals=11                 number of metals in LT chemical evolution blocks
  """

  if (len(blockd)==4):
    block=blockd
  else:
    block=blockd
    while(len(block)<4):
      block=block+" "

  if (verbose):
    print("reading block", block)
    
  blockadd=0
  blocksub=0

  if muppi==1:
    if (verbose):
      print("MUPPI format")
    blockadd=1
  if arepo==0:
    if (verbose):	
      print("Gadget format")
    blockadd=0
  if arepo==1:
    if (verbose):	
      print("Arepo format")
    blockadd=1	
  if arepo==2:
    if (verbose):
      print("Arepo extended format")
    blockadd=4	
  if no_masses==True:
    if (verbose):	
      print("No mass block present")    
    blocksub=1

 # a questo punto per MUPPI blockadd=1 blocksub=0
		 
  if parttype not in [-1,0,1,2,3,4,5]:
    raise ValueError("wrong pattern given")
  
  if os.path.exists(filename):
    curfilename = filename
  elif os.path.exists(filename+".0"):
    curfilename = filename+".0"
  else:
    raise ValueError("File not found:"+curfilename+", (initial name:"+filename+")")
  
  head = snapshot_header(curfilename)
  format = head.format

  if (verbose):
    print("SNAP FORMAT=", format, '  (2 for GADGET3)')
  swap = head.swap
  npart = head.npart
  massarr = head.massarr
  nall = head.nall
  filenum = head.filenum
  redshift = head.redshift
  time = head.time

  
  # - description of data blocks -
  # add or change blocks as needed for your Gadget version
  data_for_type = np.zeros(6,bool) # should be set to "True" below for the species for which data is stored in the data block #by doing this, the default value is False data_for_type=[False,False,False,False,False,False]
  dt = np.float32 # data type of the data in the block
  if block=="POS ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 2
  elif block=="VEL ":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 3
  elif block=="MASS":
    data_for_type[np.where(massarr==0)] = True
    block_num = 4
    dt = np.dtype(np.float32)
    if parttype>=0 and massarr[parttype]>0:   
      if (verbose):	
        print("filling masses according to massarr")   
      return np.ones(nall[parttype],dtype=dt)*massarr[parttype]
  elif block=="RHO ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 5-blocksub
  elif block=="U   ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 6-blocksub
  elif block=="NE  ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 7-blocksub
  elif block=="SFR ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 8-blocksub
  elif block=="POT ":
    data_for_type[:] = True
    dt = np.dtype(np.float32)
    block_num = 9-blocksub
##################################################
# MUPPI's blocks

  elif block=="MHOT":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 10+blockadd-blocksub
  elif block=="MCLD":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 11+blockadd-blocksub
  elif block=="MSF ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 12+blockadd-blocksub
  elif block=="EHOT":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 13+blockadd-blocksub
  elif block=="MFST":
    data_for_type[0] = True
    dt = np.dtype(np.int32)
    block_num = 14+blockadd-blocksub
  elif block=="GRAD":  # block not found
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 15+blockadd-blocksub    


##########################################

  elif block=="Zs  ":
    data_for_type[0] = True
    data_for_type[4] = True
    dt = np.dtype((np.float32,nMetals))
    block_num = 16-blocksub
  elif block=="iM  ":  #block not found
    data_for_type[4] = True
    dt = np.dtype(np.float32)
    block_num = 17-blocksub
  elif block=="AGE ":
    data_for_type[4] = True
    data_for_type[5] = True
    block_num = 18-blocksub
  elif block=="ID  ":
    data_for_type[:] = True
    dt = np.uint32
    block_num = 19-blocksub
  elif block=="BHMA":
    data_for_type[5] = True
    dt = np.dtype(np.float32)
    block_num = 20-blocksub
  elif block=="BHMD":
    data_for_type[5] = True
    dt = np.dtype(np.float32)
    block_num = 21-blocksub
  elif block=="BHPC":
    data_for_type[5] = True
    dt = np.dtype(np.int32)
    block_num = 22-blocksub
  elif block=="ACRB":
    data_for_type[5] = True
    dt = np.dtype(np.float32)
    block_num = 23-blocksub
  elif block=="GPOS":
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 24-blocksub
  elif block=="MVIR":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 25-blocksub
  elif block=="RVIR":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 26-blocksub
  elif block=="M500":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 27-blocksub
  elif block=="R500":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 28-blocksub
  elif block=="NSUB":
    data_for_type[0] = True
    dt = np.uint32
    block_num = 29-blocksub
  elif block=="MSUB":
    data_for_type[1] = True
    dt = np.dtype(np.float32)
    block_num = 30-blocksub
  elif block=="SPOS":
    data_for_type[1] = True
    dt = np.dtype((np.float32,3))
    block_num = 31-blocksub
  elif block=="SVEL":
    data_for_type[1] = True
    dt = np.dtype((np.float32,3))
    block_num = 32-blocksub
  elif block=="SCM ":
    data_for_type[1] = True
    dt = np.dtype((np.float32,3))
    block_num = 32-blocksub
  elif block=="SMST":
    data_for_type[1] = True
    dt = np.dtype((np.float32,6))
    block_num = 33-blocksub
  elif block=="RHMS":
    data_for_type[1] = True
    dt = np.dtype(np.float32)
    block_num = 34-blocksub
  elif block=="FSUB":
    data_for_type[0] = True
    dt = np.uint32
    block_num = 35-blocksub
  elif block=="GRNR":
    data_for_type[1] = True
    dt = np.uint32
    block_num = 36-blocksub
  elif block=="GACC":
    data_for_type[:] = True
    dt = np.dtype((np.float32,3))
    block_num = 37-blocksub
  elif block=="HACC":
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 38-blocksub
  elif block=="WTME":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 39-blocksub
  elif block=="WVEL":
    data_for_type[0] = True
    dt = np.dtype((np.float32,3))
    block_num = 40-blocksub
  elif block=="BHEO":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 41-blocksub
  elif block=="BHEN":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 42-blocksub
  elif block=="TEMP":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 42-blocksub
  elif block=="HOTT":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 43-blocksub
  elif block=="NOAC":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 43-blocksub
  elif block=="H2I ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 44-blocksub
  elif block=="HII ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 45-blocksub
  elif block=="HeI ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 46-blocksub
  elif block=="HeII":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 47-blocksub
  elif block=="He3 ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 48-blocksub
  elif block=="H2II":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 49-blocksub
  elif block=="HM  ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 50-blocksub
  elif block=="HD  ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 51-blocksub
  elif block=="DI  ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 51-blocksub
  elif block=="DII ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 52-blocksub
  elif block=="HeHp":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 53-blocksub
  elif block=="NH  ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 53-blocksub
  elif block=="TNAC":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 54-blocksub
  elif block=="ZsIa":
    data_for_type[0] = True
    data_for_type[4] = True
    dt = np.dtype((np.float32,nMetals))
    block_num = 55-blocksub
  elif block=="ZsAG":
    data_for_type[0] = True
    data_for_type[4] = True
    dt = np.dtype((np.float32,nMetals))
    block_num = 56-blocksub
  elif block=="MNMP":
    data_for_type[0] = True
    dt = np.dtype(np.int32)
    block_num = 57-blocksub
  elif block=="CLDX":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 58-blocksub
  elif block=="RHWP":
    data_for_type[0] = True
    dt = np.dtype(np.int32)
    block_num = 59-blocksub
  elif block=="DETI":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 60-blocksub
  elif block=="KSTP":
    data_for_type[0] = True
    dt = np.dtype(np.int32)
    block_num = 61-blocksub
  elif block=="KTME":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 62-blocksub
  elif block=="TDYN":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 63-blocksub
  elif block=="HTNH":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 64-blocksub
  elif block=="HTTH":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 65-blocksub
  elif block=="HTFH":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 66-blocksub
  elif block=="HTMH":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 67-blocksub
  elif block=="HTEH":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 68-blocksub
  elif block=="HTA ":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 69-blocksub
  elif block=="HTDT":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 70-blocksub
  elif block=="NMF ":
    data_for_type[0] = True
    dt = np.dtype(np.int32)
    block_num = 71-blocksub
  elif block=="KAGE":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 72-blocksub
  elif block=="DSTL":
    data_for_type[0] = True
    dt = np.dtype((np.float32,nMetals))
    block_num = 73-blocksub
  elif block=="DSTS":
    data_for_type[0] = True
    dt = np.dtype((np.float32,nMetals))
    block_num = 74-blocksub
  elif block=="AWND":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 75-blocksub
  elif block=="VWND":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 76-blocksub
  elif block=="AVEN":
    data_for_type[5] = True
    dt = np.dtype(np.float32)
    block_num = 77-blocksub
  elif block=="USEN":
    data_for_type[5] = True
    dt = np.dtype(np.float32)
    block_num = 78-blocksub
  elif block=="BKCK":
    data_for_type[0] = True
    dt = np.dtype(np.uint32)
    block_num = 79-blocksub
  elif block=="BHWI":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 80-blocksub
  elif block=="PID ":
    data_for_type[2] = True
    dt = np.dtype(np.uint32)
    block_num = 81-blocksub
  elif block=="SLEN":
    data_for_type[1] = True
    dt = np.dtype(np.uint32)
    block_num = 82-blocksub
  elif block=="SOFF":
    data_for_type[1] = True
    dt = np.dtype(np.uint32)
    block_num = 83-blocksub
  elif block=="MCRI":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 84-blocksub
  elif block=="RCRI":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 85-blocksub
  elif block=="M500":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 86-blocksub
  elif block=="R500":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 87-blocksub
  elif block=="M5CC":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 88-blocksub
  elif block=="R5CC":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 89-blocksub
  elif block=="M25K":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 90-blocksub
  elif block=="R25K":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 91-blocksub
  elif block=="M200":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 92-blocksub
  elif block=="R200":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 93-blocksub
  elif block=='GLEN':
    data_for_type[0] = True
    dt = np.dtype(np.uint32)
    block_num = 94-blocksub
  elif block=='GOFF':
    data_for_type[0] = True
    dt = np.dtype(np.uint32)
    block_num = 95-blocksub
  elif block=="HSML":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 96-blocksub
  elif block=="MTOT":
    data_for_type[0] = True
    dt = np.dtype(np.float32)
    block_num = 97-blocksub
  elif block=="TNGB":
    data_for_type[0] = True
    data_for_type[4] = True
    data_for_type[5] = True
    dt = np.dtype(np.int32)
    block_num = 98-blocksub
  elif block=="ACRS":
    data_for_type[4] = True
    dt = np.dtype(np.float32)
    block_num = 99-blocksub
  elif block=="RHOS":
    data_for_type[4] = True
    dt = np.dtype(np.float32)
    block_num = 100-blocksub
  elif block=="NEIS":
    data_for_type[4] = True
    dt = np.dtype(np.float32)
    block_num = 101-blocksub
  elif block=="TSTP":
    data_for_type[:] = True
    dt = np.dtype(np.float32)
    block_num = 102-blocksub
 
















# #only used for format I, when file structure is HEAD,POS,VEL,ID,ACCE
#   elif block=="ACCE":              #This is only for the PIETRONI project
#     data_for_type[:] = True        #This is only for the PIETRONI project
#     dt = np.dtype((np.float32,3))  #This is only for the PIETRONI project
#     block_num = 5                  #This is only for the PIETRONI project

#   elif block=="VOL ":
#     data_for_type[0] = True
#     block_num = 8-blocksub 
#   elif block=="CMCE":
#     data_for_type[0] = True
#     dt = np.dtype((np.float32,3))
#     block_num = 9-blocksub 
#   elif block=="AREA":
#     data_for_type[0] = True
#     block_num = 10-blocksub
#   elif block=="NFAC":
#     data_for_type[0] = True
#     dt = np.dtype(np.int64)        #depends on code version, most recent hast int32, old MyIDType	
#     block_num = 11-blocksub
#   elif block=="NE  ":
#     data_for_type[0] = True
#     block_num = 8+blockadd-blocksub
#   elif block=="NH  ":
#     data_for_type[0] = True
#     block_num = 9+blockadd-blocksub
#   elif block=="HSML":
#     data_for_type[0] = True
#     block_num = 10+blockadd-blocksub
#   elif block=="Z   ": # aggiunta
#     data_for_type[0] = True
#     data_for_type[4] = True
#     block_num = 13+blockadd-blocksub

  else:
    raise ValueError("Sorry! Block type"+ str(block)+ "not known!")
  # - end of block description -
    
  actual_data_for_type = np.copy(data_for_type)  
  if parttype >= 0:
    actual_data_for_type[:] = False
    actual_data_for_type[parttype] = True
    if data_for_type[parttype]==False:
      raise ValueError("Error: no data for specified particle type"+str( parttype)+ "in the block"+str( block   ))

  elif block=="MASS":
    actual_data_for_type[:] = True  
    
  allpartnum = np.int64(0)
  species_offset = np.zeros(6,np.int64)
  for j in range(6):
    species_offset[j] = allpartnum
    if actual_data_for_type[j]:
      allpartnum += nall[j]
    
  for i in range(filenum): # main loop over files
    if filenum>1:
      curfilename = filename+"."+str(i)
      
    if i>0:
      head = snapshot_header(curfilename)
      npart = head.npart  
      del head

    if block=='BHMA' and  npart[5]==0:
      data = np.empty(allpartnum,dt)
      continue
      
    curpartnum = np.int32(0)
    cur_species_offset = np.zeros(6,np.int64)
    for j in range(6):
      cur_species_offset[j] = curpartnum
      if data_for_type[j]:
        curpartnum += npart[j]
    
    if parttype>=0:
      actual_curpartnum = npart[parttype]      
      add_offset = cur_species_offset[parttype] 
    else:
      actual_curpartnum = curpartnum
      add_offset = np.int32(0)
      
    offset,blocksize = find_block(curfilename,format,swap,block,block_num)

    
    if i==0: # fix data type for ID if long IDs are used
      if block=="ID  ":
        if blocksize == np.dtype(dt).itemsize*curpartnum * 2:
          dt = np.uint64 


    if np.dtype(dt).itemsize*curpartnum != blocksize:
       print(" ")
       print("--------")
       print(curpartnum, np.dtype(dt).itemsize)
       raise ValueError("something wrong with blocksize! expected ="+str(np.dtype(dt).itemsize*curpartnum)+" actual ="+str(blocksize))

    
    f = open(curfilename,'rb')
    f.seek(offset + add_offset*np.dtype(dt).itemsize, os.SEEK_CUR)  
    curdat = np.fromfile(f,dtype=dt,count=actual_curpartnum) # read data
    f.close()  
    if swap:
      curdat.byteswap(True)  
      
    if i==0:
      data = np.empty(allpartnum,dt)
    
    for j in range(6):
      if actual_data_for_type[j]:
        if block=="MASS" and massarr[j]>0: # add mass block for particles for which the mass is specified in the snapshot header
          data[species_offset[j]:species_offset[j]+npart[j]] = massarr[j]
        else:
          if parttype>=0:
            data[species_offset[j]:species_offset[j]+npart[j]] = curdat
          else:
            data[species_offset[j]:species_offset[j]+npart[j]] = curdat[cur_species_offset[j]:cur_species_offset[j]+npart[j]]
        species_offset[j] += npart[j]

    del curdat

  if physical_velocities and block=="VEL " and redshift!=0:
    data *= math.sqrt(time)

  return data
  
# ----- list all data blocks in a format 2 snapshot file -----

def list_format2_blocks(filename):
  if os.path.exists(filename):
    curfilename = filename
  elif os.path.exists(filename+".0"):
    curfilename = filename+".0"
  else:
    raise ValueError("File not found:"+filename)
  
  head = snapshot_header(curfilename)
  format = head.format
  swap = head.swap
  del head
  
  print('GADGET FORMAT ',format)
  if (format != 2):
    print("#   OFFSET   SIZE")
  else:            
    print("#   BLOCK   OFFSET   SIZE")
  print("-------------------------")
  
  find_block(curfilename, format, swap, "XXXX", 0, only_list_blocks=True)
  
  print("-------------------------")

def read_gadget_header(filename):
  if os.path.exists(filename):
    curfilename = filename
  elif os.path.exists(filename+".0"):
    curfilename = filename+".0"
  else:
    raise ValueError("File not found:"+filename)

  head=snapshot_header(curfilename)
  print('npar=',head.npart)
  print('nall=',head.nall)
  print('a=',head.time)
  print('z=',head.redshift)
  print('masses=',head.massarr*1e10,'Msun/h')
  print('boxsize=',head.boxsize,'kpc/h')
  print('filenum=',head.filenum)
  print('cooling=',head.cooling)
  print('Omega_m,Omega_l=',head.omega_m,head.omega_l)
  print('h=',head.hubble,'\n')
  
  rhocrit=2.77536627e11 #h**2 M_sun/Mpc**3
  rhocrit=rhocrit/1e9 #h**2M_sun/kpc**3
  
  Omega_DM=head.nall[1]*head.massarr[1]*1e10/(head.boxsize**3*rhocrit)
  print('DM mass=',head.massarr[1]*1e10,'Omega_DM=',Omega_DM)
  if head.nall[2]>0 and head.massarr[2]>0:
    Omega_NU=head.nall[2]*head.massarr[2]*1e10/(head.boxsize**3*rhocrit)
    print('NU mass=',head.massarr[2]*1e10,'Omega_NU=',Omega_NU)
    print('Sum of neutrino masses=',Omega_NU*head.hubble**2*94.1745,'eV')


def detect_block(filename, blockd):
  """
      USE: detect_block('snapfile','BLC")
      
      returns True if the block is present, False otherwise
      auto-recognizes format and swap
  """
  if (not os.path.exists(filename)):
      raise ValueError("File not found:"+filename)


  f = open(filename,'rb')    
  blocksize = np.fromfile(f,dtype=np.int32,count=1)
  if blocksize[0] == 8:
    swap = 0
    format = 2
  elif blocksize[0] == 256:
    swap = 0
    format = 1  
  else:
    blocksize.byteswap(True)
    if blocksize[0] == 8:
      swap = 1
      format = 2
    elif blocksize[0] == 256:
      swap = 1
      format = 1
    else:
      raise ValueError("incorrect file format encountered when reading header of"+ filename)
  f.close()

  if (len(blockd)==4):
    block=blockd
  else:
    block=blockd
    while(len(block)<4):
      block=block+" "

  f = open(filename,'rb')
  f.seek(0, os.SEEK_END)
  filesize = f.tell()
  f.seek(0, os.SEEK_SET)
  
  found = False
  curblock_num = 1
  while ((not found) and (f.tell()<filesize)):
    if format==2:
      f.seek(4, os.SEEK_CUR)
      curblock = f.read(4)
#      print curblock_num, block, curblock
      if (block == curblock):
        found = True
      f.seek(8, os.SEEK_CUR)  
    else:
      if curblock_num==block_num:
        found = True

    curblocksize = (np.fromfile(f,dtype=np.uint32,count=1))[0]
    if swap:
      curblocksize = curblocksize.byteswap()

#    print " ------ ",curblocksize
    if found:
      blocksize = curblocksize
      offset = f.tell()
      f.close()
      return True
    else:
      f.seek(curblocksize, os.SEEK_CUR)
      blocksize_check = (np.fromfile(f,dtype=np.uint32,count=1))[0]
      if swap: blocksize_check = blocksize_check.byteswap()
      if (curblocksize != blocksize_check):
        raise ValueError("something wrong")
      curblock_num += 1

  f.close()
  return False
