c---------------------------------------------------------------------
c                         PROGRAM RENUM2          ECMeng   2/98
c  --interactive
c  --renumbers residues sequentially, removing insertion codes and chain ID's
c  --allows user specification of starting residue number
c  --does not renumber the atoms
c  --does not save TER cards (commented out)
c---------------------------------------------------------------------
c
      character*80 pdblin, pdbfil, pdbout
      character*5 resid1, resid2
      integer count
c
      write (6, *) 'enter name of pdb input file:'
      read (5, 1000) pdbfil
 1000 format (A66)
      write (6, 1000) pdbfil
      write (6, *) 'enter name of output file:'
      read (5, 1000) pdbout
      write (6, 1000) pdbout
      write (6, *) 'number to start with?'
      read (5,*) count
      open (unit=1, file=pdbfil, status='old')
      open (unit=2, file=pdbout, status='new')
c
c     count = 1
   50 read (1, 1000, end=500) pdblin
      if (pdblin(1:4).eq.'ATOM'.or.pdblin(1:4).eq.'HETA') then
        read (pdblin, '(22x,a5)') resid1
        write (2,1001) pdblin(1:20),count,pdblin(28:66)
 1001   format (a20,2x,i4,a39)
      else if (pdblin(1:3).eq.'TER') then
        write (2,1002) 'TER'
        count = count + 1
 1002   format (a3)
      endif
      goto 50
  500 continue
      write (6, *) 'number of last residue written = ', count
      close (1)
      close (2)
      end
c---------------------------------------------------------------------

