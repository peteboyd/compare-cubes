      program create_esp
      implicit none
!-----General parameters
      integer n_dim, n_atoms_max, n_atoms_type_max 
       parameter(n_dim=3,n_atoms_max=5000,n_atoms_type_max=171)

      integer n_atoms1,ngx1,ngy1,ngz1,
     & n_atoms2,ngx2,ngy2,ngz2,
     & n_atoms,ngx,ngy,ngz,
     & ngpt
      double precision axis_zero(n_dim), axis_vector(n_dim,n_dim),
     & real_box_vector(n_dim,n_dim), 
     & Box_volume, recip_box_vector(n_dim,n_dim)
      common/dbl_general/axis_zero, axis_vector,
     & real_box_vector, Box_volume, recip_box_vector

!-----Mathematical constants
      double precision pi, TWOPI
      parameter(pi=3.14159265,TWOPI=6.2831853)
      double precision vdw_fact
      common/dbl_input/vdw_fact
      double precision, allocatable :: vdw_radii(:)
!-----Output results in a cube file?
      character*1, allocatable :: rootname(:)
      character(len=70) file_name
      common/char_input/file_name

c     WARNING: Setting lenrec too high (eg 256) results in
c     undefined behaviour when performing string manipulation
c     in the cif file.
      integer, parameter :: lenrec = 120
      integer, parameter :: ncube = 205
      real(8), parameter :: angs2bohr = 1.889725989d0
     
      integer, allocatable :: atom_number(:),atom_index(:)
      real(8), allocatable :: eps(:),sig(:)
      character*2, allocatable :: atom(:)
      character*3, allocatable :: atmname(:)
      character*1 title(lenrec)
      integer i_atom, seed, i_dim, i, j, k, i_neigh_k, n_loop,
     & i_neigh_j, i_neigh_i, i_extra, i_loop, j_dim, loopcount
    
      integer idum, ind 
      double precision, allocatable, dimension(:,:) :: atom_pos
      double precision, allocatable, dimension(:,:) :: atom_pos_frac
      real(8) alen,blen,clen,alph,beta,gamma
      character*1 record(lenrec)
      double precision, allocatable, dimension(:,:,:) :: V_pot1,V_pot2,
     &V_diff
      integer, allocatable, dimension(:,:,:) :: V_flag
      real(8) time,ftime,etime
    
      double precision delta_dist(n_dim), delta_fdist(n_dim),dist,
     &grid_pos_frac(n_dim),DO_FACTOR
      parameter(DO_FACTOR=1.0)
      character(len=100) :: cube1,cube2
      save record,atom,atmname,rootname
      save atom_number,atom_index

      double precision Phi_sum1,Phi_sum2,diff,diff2,MAD,SSE 
      
      if (iargc().ne.2) then
          write(*,*)"Usage: compare_cubes [cube file 1] [cube file 2]"
          call exit(1)
      end if
c     start timing
      call timchk(0, time)
c     get filenames from command line 
      call getarg(1, cube1)
      call getarg(2, cube2)

      call scancube(ncube,cube1,n_atoms1,ngx1,ngy1,ngz1)
      call scancube(ncube,cube2,n_atoms2,ngx2,ngy2,ngz2)
      if(n_atoms1.ne.n_atoms2)then
          write(*,*)"ERROR - the cube file atom counts don't match!"
          stop
      else if((ngx1.ne.ngx2).or.(ngy1.ne.ngy2).or.(ngz1.ne.ngz2))then
          write(*,*)"ERROR - the cube files have different number of 
     & grid points!"
          stop
      end if

      n_atoms = n_atoms1
      ngx = ngx1
      ngy = ngy1
      ngz = ngz1
      allocate(atom_index(n_atoms))
      allocate(atom_number(n_atoms))
      allocate(atom_pos(n_atoms,n_dim))
      allocate(atom_pos_frac(n_atoms,n_dim))
      allocate(vdw_radii(n_atoms))
      allocate(V_pot1(ngx,ngy,ngz))
      allocate(V_pot2(ngx,ngy,ngz))
      allocate(V_diff(ngx,ngy,ngz))
      allocate(V_flag(ngx,ngy,ngz))
      call scancubeatoms(ncube,cube1,n_atoms1)
!     compute fractional positions of atoms
      call process_box(n_atoms,ngx,ngy,ngz)
      call VDW_radii_array(vdw_radii,n_atoms)
      ! init to 0
      do i=1,ngx
        do j=1,ngy
          do k=1,ngz
            V_pot1(i,j,k)=0.d0
            V_pot2(i,j,k)=0.d0
            V_diff(i,j,k)=0.d0
            V_flag(i,j,k)=1
            if(DO_FACTOR.ne.0)then
              ! flag if not to compute ESP at this grid point
              grid_pos_frac(1) = dble(i-1)/dble(ngx)
              grid_pos_frac(2) = dble(j-1)/dble(ngy)
              grid_pos_frac(3) = dble(k-1)/dble(ngz)
              do i_atom=1, n_atoms
                delta_fdist(1)=grid_pos_frac(1)-atom_pos_frac(i_atom,1)
                delta_fdist(2)=grid_pos_frac(2)-atom_pos_frac(i_atom,2)
                delta_fdist(3)=grid_pos_frac(3)-atom_pos_frac(i_atom,3)
                ! shift by pbc
                delta_fdist(1)=delta_fdist(1)-dble(nint(delta_fdist(1)))
                delta_fdist(2)=delta_fdist(2)-dble(nint(delta_fdist(2)))
                delta_fdist(3)=delta_fdist(3)-dble(nint(delta_fdist(3)))

                ! compute cartesian
                dist=0.d0
                do i_dim=1, 3
                  delta_dist(i_dim) =
     &              delta_fdist(1)*real_box_vector(1,i_dim)+
     &              delta_fdist(2)*real_box_vector(2,i_dim)+
     &              delta_fdist(3)*real_box_vector(3,i_dim)
                !  atom_pos_tmp(i_dim) = atom_pos(i_atom,i_dim)
                !  delta_dist(i_dim) = grid_pos(i_dim) - atom_pos_tmp(i_dim)
                  dist = dist + delta_dist(i_dim)**2
                end do
                ! check for nearby atoms
                dist = sqrt(dist)
                if(dist.le.(vdw_radii(i_atom)+DO_FACTOR))then
                  V_flag(i,j,k)=0
                endif
              end do
            endif

          enddo
        enddo
      enddo
      call readcube(ncube,cube1,n_atoms,V_pot1,ngx,ngy,ngz)
      call readcube(ncube,cube2,n_atoms,V_pot2,ngx,ngy,ngz)
      
      call getroot(cube1,ind)
      call timchk(0,etime)

!     comupte the difference squared between the cube grids.
      Phi_sum1 = 0.d0
      Phi_sum2 = 0.d0
      MAD = 0.d0
      SSE = 0.d0
      ngpt = 0
      do i=1, ngx
       do j=1, ngy
        do k=1, ngz
          if(V_flag(i,j,k).eq.1)then
            Phi_sum1 = Phi_sum1 + V_pot1(i,j,k) 
            Phi_sum2 = Phi_sum2 + V_pot2(i,j,k)
            diff = V_pot1(i,j,k) - V_pot2(i,j,k)
            diff2 = diff*diff
            V_diff(i,j,k) = diff
            SSE = SSE+diff2
            MAD = MAD + diff
            ngpt = ngpt + 1
          endif
        end do 
       end do
      end do

      call timchk(1,etime)

      call timchk(0,ftime)
      call writecube(ind,V_diff,n_atoms,ngx,ngy,ngz)

      call timchk(1,ftime)
      call timchk(1,time)

      write(*,'(a,f20.5)')"Total Sum Sq. Error :",SSE
      write(*,'(a,f20.5)')"Total Mean Abs. Dev.:",MAD/dble(ngpt)
      write(*,'(a,f20.5,a)')"ESP comparison time :",etime," seconds."
      write(*,'(a,f20.5,a)')"File writing time   :",ftime," seconds."
      write(*,'(a,f20.5,a)')"Total wall time     :",time," seconds."
      contains 
      
      subroutine writecube(ind,V_pot,n_atoms,ngx,ngy,ngz)
c***********************************************************************
c     
c     routine to write a cube file from the potential generated
c     in the main program
c
c***********************************************************************
      implicit none
      integer i_dim, i_atoms, ind, n_loop, i_extra
     & i, j, k, n_atoms, i_atom, i_loop,
     & ngx, ngy, ngz
      double precision V_pot(ngx,ngy,ngz)
      character(len=ind+10) cubefile
      character(len=ind) root
      character*4 p
      do i=1, ind
        root(i:i) = rootname(i)
      end do
      cubefile=root//"_diff.cube"

      open(10,file=cubefile,status="new")
!----Write cube file
      write(10,*)"Cube file generated with compare-cubes"
      write(10,*)"*****************************************"

!-----Writing the number of atoms and origin
      write(10,'(i5,3f12.6)') n_atoms,
     &      (axis_zero(i_dim), i_dim=1,n_dim)

!-----Writing the voxels arrray
      write(10,'(i5,3f12.6)')ngx, 
     &      (axis_vector(1,i_dim), i_dim=1,3)
      write(10,'(i5,3f12.6)')ngy, 
     &      (axis_vector(2,i_dim), i_dim=1,3)
      write(10,'(i5,3f12.6)')ngz, 
     &      (axis_vector(3,i_dim), i_dim=1,3)

      do i_atom=1, n_atoms
c        write(*,*)atom_number(i_atom),atom_index(i_atom),
c     &atom_pos(i_atom,1:3)
        write(10,'(2i5,3f12.6)')atom_number(i_atom), atom_index(i_atom),
     & (atom_pos(i_atom,i_dim), i_dim=1, 3) 
      end do
      n_loop = floor(real(ngz)/6.d0)
      i_extra = mod(ngz,6)
      write(p,'(i1)') i_extra
c      write(*,*) p, i_extra, n_loop

      do i=1, ngx
       do j=1, ngy
        if(i_extra.ne.0) then    
         do i_loop=1, n_loop 
          write(10,'(6e13.5)') (V_pot(i,j,k+(i_loop-1)*6), k=1, 6) 
         end do          
         write(10,'('//trim(p)//'e13.5)') (V_pot(i,j,k+(i_loop-1)*6), 
     &   k=1, i_extra)         
        else
         do i_loop=1, n_loop 
          write(10,'(6e13.5)') (V_pot(i,j,k+(i_loop-1)*6), k=1, 6)
         end do
        end if          
       end do
      end do
      close(10)
      end subroutine writecube

      real(8) function atmmass(i) 

c***********************************************************************
c     
c     function to generate an atom label from an atomic 
c     number 
c
c***********************************************************************
      implicit none
      character*2 i
      if (i.eq.'H ')then
        atmmass=1.00794
      elseif(i.eq.'He')then
        atmmass=4.002602
      elseif(i.eq.'Li')then
        atmmass=6.941
      elseif(i.eq.'Be')then
        atmmass=9.012182
      elseif(i.eq.'B ')then
        atmmass=10.811
      elseif(i.eq.'C ')then
        atmmass=12.0107
      elseif(i.eq.'N ')then
        atmmass=14.00674
      elseif(i.eq.'O ')then
        atmmass=15.9994
      elseif(i.eq.'F ')then
        atmmass=18.9984032
      elseif(i.eq.'Ne')then
        atmmass=20.1797
      elseif(i.eq.'Na')then
        atmmass=22.989770
      elseif(i.eq.'Mg')then
        atmmass=24.3050
      elseif(i.eq.'Al')then
        atmmass=26.981538
      elseif(i.eq.'Si')then
        atmmass=28.0855
      elseif(i.eq.'P ')then
        atmmass=30.973761
      elseif(i.eq.'S ')then
        atmmass=32.066
      elseif(i.eq.'Cl')then
        atmmass=35.4527
      elseif(i.eq.'Ar')then
        atmmass=39.948
      elseif(i.eq.'K ')then
        atmmass=39.0983
      elseif(i.eq.'Ca')then
        atmmass=40.078
      elseif(i.eq.'Sc')then
        atmmass=44.955910
      elseif(i.eq.'Ti')then
        atmmass=47.867
      elseif(i.eq.'V ')then
        atmmass=50.9415
      elseif(i.eq.'Cr')then
        atmmass=51.9961
      elseif(i.eq.'Mn')then
        atmmass=54.938049
      elseif(i.eq.'Fe')then
        atmmass=55.845
      elseif(i.eq.'Co')then
        atmmass=58.9332
      elseif(i.eq.'Ni')then
        atmmass=58.6934
      elseif(i.eq.'Cu')then
        atmmass=63.546
      elseif(i.eq.'Zn')then
        atmmass=65.39
      elseif(i.eq.'Ga')then
        atmmass=69.723
      elseif(i.eq.'Ge')then
        atmmass=72.61
      elseif(i.eq.'As')then
        atmmass=74.9216
      elseif(i.eq.'Se')then
        atmmass=78.96
      elseif(i.eq.'Br')then
        atmmass=79.904
      elseif(i.eq.'Kr')then
        atmmass=83.80
      elseif(i.eq.'Y ')then
        atmmass=88.90585
      elseif(i.eq.'Zr')then
        atmmass=91.224
      elseif(i.eq.'Nb')then
        atmmass=92.90638
      elseif(i.eq.'Mo')then
        atmmass=95.94
      elseif(i.eq.'Ru')then
        atmmass=101.07
      elseif(i.eq.'Rh')then
        atmmass=102.90550
      elseif(i.eq.'Pd')then
        atmmass=106.42
      elseif(i.eq.'Ag')then
        atmmass=107.8682
      elseif(i.eq.'Cd')then
        atmmass=112.411
      elseif(i.eq.'In')then
        atmmass=114.818
      elseif(i.eq.'Sn')then
        atmmass=118.710
      elseif(i.eq.'Sb')then
        atmmass=121.760
      elseif(i.eq.'Te')then
        atmmass=127.760
      elseif(i.eq.'I ')then
        atmmass=126.90447
      elseif(i.eq.'Xe')then
        atmmass=131.29
      elseif(i.eq.'Ba')then
        atmmass=137.327
      elseif(i.eq.'Hf')then
        atmmass=178.49
      elseif(i.eq.'Ta')then
        atmmass=180.9479
      elseif(i.eq.'W ')then
        atmmass=183.84
      elseif(i.eq.'Re')then
        atmmass=186.207
      elseif(i.eq.'Os')then
        atmmass=190.23
      elseif(i.eq.'Ir')then
        atmmass=192.217
      elseif(i.eq.'Pt')then
        atmmass=195.078
      elseif(i.eq.'Au')then
        atmmass=196.96655
      elseif(i.eq.'Hg')then
        atmmass=200.59
      elseif(i.eq.'Tl')then
        atmmass=204.3833
      elseif(i.eq.'Pb')then
        atmmass=207.2
      endif
      return
      end function atmmass

      integer function intstr(word,len,lst)

c***********************************************************************
c     
c     function for extracting integers from a 
c     character string
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c
c***********************************************************************
      
      implicit none

      logical flag,count,final
      character*1 n,word,ksn
      integer lst,len,j,isn

      dimension n(0:9),word(len)
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      lst=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      count=.false.
      
      do while(lst.lt.len.and.(.not.final))

        lst=lst+1
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            count=.true.
            flag=.true.
            
          endif
          
        enddo

        if(count.and.(.not.flag))final=.true.
        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

      enddo

      intstr=isn*intstr

      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      return
      end function intstr

      real(8) function dblstr(word,len,lst)

c***********************************************************************
c     
c     Function for extracting double precisions from a 
c     character string. 
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c     
c***********************************************************************
      
      implicit none
      
      character*1 n,word,ksn,dot,d,e
      logical flag,ldot,start,final
      integer len,lst,iexp,idum,i,j,fail
      real(8) sn,ten,one
      dimension n(0:9),word(len)
      character*1, allocatable :: work(:)

      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      allocate(work(len),stat=fail)

      lst=0
      sn=1.d0
      ksn='+'
      ten=10.d0
      one=1.d0
      
      dblstr=0.d0
      iexp=0
      idum=0
      start=.false.
      ldot=.false.
      final=.false.

      do while(lst.lt.len.and.(.not.final))
        
        lst=lst+1
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr=ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst))then
          
          flag=.true.
          ten=1.d0
          ldot=.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        if(ldot) one=one/10.d0
        ksn=word(lst)
        if(ksn.eq."D")ksn="d"
        if(ksn.eq."E")ksn="e"
        
        if(start)then
          
          if(d.eq.ksn.or.e.eq.ksn)then
            
            do i=1,len-lst
              work(i)=word(i+lst)
            enddo
            iexp=intstr(work,len-lst,idum)
            final=.true.

          endif

          if(.not.flag)final=.true.
          
        endif
        
      enddo
      
      dblstr=sn*dblstr*(10.d0**iexp)
      lst=lst+idum
      
      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      deallocate(work,stat=idum)

      return
      end function dblstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      subroutine scancube(ncube,cubefile,n_atoms,ngx,ngy,ngz)
c*********************************************************************
c
c     routine to scan cube file for array allocation 
c
c*********************************************************************
      implicit none
      integer ncube,n_atoms,ngx,ngy,ngz
      character*100 cubefile

      open(ncube,file=cubefile,status="old")
!-----Reading the first two comments
      read(ncube,*)
      read(ncube,*)
!-----Reading the number of atoms and origin
      read(ncube,*) n_atoms, (axis_zero(i_dim), i_dim=1,3)
!-----Reading the voxels arrray
      read(ncube,*) ngx, (axis_vector(1,i_dim), i_dim=1,3)
      read(ncube,*) ngy, (axis_vector(2,i_dim), i_dim=1,3)
      read(ncube,*) ngz, (axis_vector(3,i_dim), i_dim=1,3)
!-----Reading atomic info
      close(ncube)

      return
      end subroutine scancube
      
      subroutine scancubeatoms(ncube,cubefile,n_atoms)
c*********************************************************************
c
c     routine to scan cube file for atomic info
c
c*********************************************************************
      implicit none
      integer ncube,n_atoms,i_atom
      character*100 cubefile,junk

      open(ncube,file=cubefile,status="old")
!-----Reading the first two comments
      read(ncube,*)
      read(ncube,*)
!-----Reading the number of atoms and origin
      read(ncube,*) 
!-----Reading the voxels arrray
      read(ncube,*) 
      read(ncube,*) 
      read(ncube,*) 
!-----Reading atomic info
      do i_atom=1, n_atoms
       read(ncube,*) atom_number(i_atom), atom_index(i_atom),
     & (atom_pos(i_atom,i_dim), i_dim=1, 3)
      end do
      close(ncube)

      return
      end subroutine scancubeatoms

      subroutine readcube(ncube,cubefile,n_atoms,V_pot,ngx,ngy,ngz)
c*********************************************************************
c
c     routine to read cif files and store all the values in memory 
c
c*********************************************************************
      implicit none
      integer ncube,n_atoms,n_loop,i_extra
      integer ngx,ngy,ngz,i,j,k
      real(8) xcoord,ycoord,zcoord
      double precision V_pot(ngx,ngy,ngz)
      character*8 atm,atmn,junk
      character*100 cubefile
      character*4 p

      open(ncube,file=cubefile,status="old")
!-----Reading the first two comments
      read(ncube,*)
      read(ncube,*)
!-----Reading the number of atoms and origin
      read(ncube,*) 
!-----Reading the voxels arrray
      read(ncube,*) 
      read(ncube,*)
      read(ncube,*)
!-----Reading atomic info
      do i_atom=1, n_atoms
       read(ncube,*) 
      end do
!-----Creating the format for reading the ESP field
      n_loop = floor(real(ngz)/6)
      i_extra = mod(ngz,6)
      write(p,'(i1)') i_extra
!-----Reading voxels info
      do i=1, ngx
       do j=1, ngy
        if(i_extra.ne.0) then
         do i_loop=1, n_loop
          read(ncube,'(6e13.5)') (V_pot(i,j,k+(i_loop-1)*6), k=1, 6)
         end do
         read(ncube,'('//trim(p)//'e13.5)') (V_pot(i,j,k+(i_loop-1)*6),
     &   k=1, i_extra)
        else
         do i_loop=1, n_loop
          read(ncube,'(6e13.5)') (V_pot(i,j,k+(i_loop-1)*6), k=1, 6)
         end do
        end if
       end do
      end do
      
      close(ncube) 
      return
      end subroutine readcube

      logical function findstring(seek,string,here)

c***********************************************************************
c     
c     routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c     
c***********************************************************************

      implicit none

      integer i,n,m,here
      character*(*) seek
      character*1 string(lenrec)

      m=lenrec
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end function findstring
      
      subroutine copystring(oldstr,newstr,length)

c***********************************************************************
c     
c     routine to copy one string into another
c     
c***********************************************************************

      implicit none

      character*1 newstr(*),oldstr(*)
      integer i,length

      do i=1,length

        newstr(i)=oldstr(i)

      enddo

      return
      end subroutine copystring
      
      subroutine getword(word,string,len1,len2)

c***********************************************************************
c     
c     routine to fetch an 8 character word from a string
c     while ignoring leading blanks
c     
c***********************************************************************

      implicit none

      logical final
      character*8 word
      integer len1,len2,i,j,k
      character*1 wrdseq(len1),string(len2)
c      character*1 word(len1),string(len2)
      do i=1,len1
        wrdseq(i)=' '
c         word(i)=' '
      enddo

      i=0
      k=0
      final=.false.
      
      do while(.not.final.and.i.lt.len2)
        
        i=i+1
        
        if(string(1).eq.' ')then
          
          if(k.gt.0)final=.true.
          
        else
          
          k=k+1
          wrdseq(k)=string(1)
c          word(k)=string(1)
          if(k.eq.len1)final=.true.

        endif
        
        do j=1,len2-1
          
          string(j)=string(j+1)
          
        enddo
        
        string(len2)=' '
          
      enddo
      
      word=mkwd8(wrdseq)

      return
      end subroutine getword
      
      subroutine getrec(safe,ifile)

c*********************************************************************
c     
c      subroutine to read a character string on one node
c     
c*********************************************************************

      implicit none
      
      logical safe

      character*150 line
      integer ifile,i
      
      safe=.true.
      
      read(ifile,'(a150)',end=100)line
 
      do i=1,lenrec

        record(i)=line(i:i)
         
      enddo
             
      return
        
  100 safe=.false.
                  
      end subroutine getrec


      subroutine lowcase(string,length)

c***********************************************************************
c     
c     routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c     
c***********************************************************************

      implicit none

      character*1 string(*)
      character*1 letter
      integer i,length

      do i=1,min(255,length)

        letter=string(i)

        if(letter.eq.'A')then
          letter='a'
        else if(letter.eq.'B')then
          letter='b'
        else if(letter.eq.'C')then
          letter='c'
        else if(letter.eq.'D')then
          letter='d'
        else if(letter.eq.'E')then
          letter='e'
        else if(letter.eq.'F')then
          letter='f'
        else if(letter.eq.'G')then
          letter='g'
        else if(letter.eq.'H')then
          letter='h'
        else if(letter.eq.'I')then
          letter='i'
        else if(letter.eq.'J')then
          letter='j'
        else if(letter.eq.'K')then
          letter='k'
        else if(letter.eq.'L')then
          letter='l'
        else if(letter.eq.'M')then
          letter='m'
        else if(letter.eq.'N')then
          letter='n'
        else if(letter.eq.'O')then
          letter='o'
        else if(letter.eq.'P')then
          letter='p'
        else if(letter.eq.'Q')then
          letter='q'
        else if(letter.eq.'R')then
          letter='r'
        else if(letter.eq.'S')then
          letter='s'
        else if(letter.eq.'T')then
          letter='t'
        else if(letter.eq.'U')then
          letter='u'
        else if(letter.eq.'V')then
          letter='v'
        else if(letter.eq.'W')then
          letter='w'
        else if(letter.eq.'X')then
          letter='x'
        else if(letter.eq.'Y')then
          letter='y'
        else if(letter.eq.'Z')then
          letter='z'
        endif

        string(i)=letter

      enddo

      return
      end subroutine lowcase

      subroutine strip(string,imax)

c***********************************************************************
c     
c     Routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c***********************************************************************

      implicit none

      integer i,imax,j
      character*1 string(imax)
      do i=1,imax
    
        if(string(1).eq.' ')then

          do j=1,imax-1

            string(j)=string(j+1)

          enddo

          string(imax)=' '

        endif

      enddo

      return
      end subroutine strip
      
      character*8 function mkwd8(string)

c***********************************************************************
c     
c     Routine to make an 8 character word from a string
c
c***********************************************************************

      implicit none

      integer i
      character*1 string(*)
      
      do i=1,8
         mkwd8(i:i)=string(i)
      enddo
      
      return
      end function mkwd8
      
      integer function atmnumber(i) 
c*******************************************************************
c     generate an atomic number from an atomic mass 
c     EDIT (pb 09/01/13): this function reads the mass reported
c     on a standard periodic table and assigns an atomic number.
c     You will run into problems if you are using atomic masses
c     of isotopes in the FIELD file.
c*******************************************************************
      implicit none
      real(8) i
      if ((i.ge.0.0).and.(i.le.1.5))then
        atmnumber=1
      elseif((i.ge.3.9).and.(i.le.4.5))then
        atmnumber=2
      elseif((i.ge.6.5).and.(i.le.7.1))then
        atmnumber=3
      elseif((i.ge.8.9).and.(i.le.9.5))then
        atmnumber=4
      elseif((i.ge.10.5).and.(i.le.11.1))then
        atmnumber=5
      elseif((i.ge.11.9).and.(i.le.12.5))then
        atmnumber=6
      elseif((i.ge.13.9).and.(i.le.14.5))then
        atmnumber=7
      elseif((i.ge.15.5).and.(i.le.16.1))then
        atmnumber=8
      elseif((i.ge.18.5).and.(i.le.19.1))then
        atmnumber=9
      elseif((i.ge.19.9).and.(i.le.20.5))then
        atmnumber=10
      elseif((i.ge.22.5).and.(i.le.23.1))then
        atmnumber=11
      elseif((i.ge.23.9).and.(i.le.24.5))then
        atmnumber=12
      elseif((i.ge.26.5).and.(i.le.27.1))then
        atmnumber=13
      elseif((i.ge.27.9).and.(i.le.28.5))then
        atmnumber=14
      elseif((i.ge.30.5).and.(i.le.31.1))then
        atmnumber=15
      elseif((i.ge.31.9).and.(i.le.32.5))then
        atmnumber=16
      elseif((i.ge.34.9).and.(i.le.36.1))then
        atmnumber=17
c     Ar (18) has mass range that overlaps with Ca (20). Be careful 
c     with mass rounding here!
      elseif((i.ge.39.5).and.(i.le.39.9999))then
        atmnumber=18
      elseif((i.ge.38.9).and.(i.le.39.4))then
        atmnumber=19
      elseif((i.ge.40.0).and.(i.le.40.5))then
        atmnumber=20
      elseif((i.ge.44.5).and.(i.le.45.1))then
        atmnumber=21
      elseif((i.ge.47.5).and.(i.le.48.1))then
        atmnumber=22
      elseif((i.ge.50.5).and.(i.le.51.1))then
        atmnumber=23
      elseif((i.ge.51.5).and.(i.le.52.1))then
        atmnumber=24
      elseif((i.ge.54.5).and.(i.le.55.1))then
        atmnumber=25
      elseif((i.ge.55.5).and.(i.le.56.1))then
        atmnumber=26
c     Co (27) and Ni (28) have very close mass ranges
      elseif((i.ge.58.76).and.(i.le.59.1))then
        atmnumber=27
      elseif((i.ge.58.5).and.(i.le.59.75))then
        atmnumber=28
      elseif((i.ge.62.9).and.(i.le.64.1))then
        atmnumber=29
      elseif((i.ge.64.9).and.(i.le.66.1))then
        atmnumber=30
      elseif((i.ge.69.5).and.(i.le.70.1))then
        atmnumber=31
      elseif((i.ge.72.5).and.(i.le.73.1))then
        atmnumber=32
      elseif((i.ge.74.5).and.(i.le.75.1))then
        atmnumber=33
      elseif((i.ge.78.5).and.(i.le.79.1))then
        atmnumber=34
      elseif((i.ge.79.5).and.(i.le.80.1))then
        atmnumber=35
      elseif((i.ge.83.5).and.(i.le.84.1))then
        atmnumber=36
      elseif((i.ge.84.9).and.(i.le.86.1))then
        atmnumber=37
      elseif((i.ge.87.5).and.(i.le.88.1))then
        atmnumber=38
      elseif((i.ge.88.5).and.(i.le.89.1))then
        atmnumber=39
      elseif((i.ge.90.9).and.(i.le.91.5))then
        atmnumber=40
      elseif((i.ge.92.5).and.(i.le.93.1))then
        atmnumber=41
      elseif((i.ge.95.5).and.(i.le.96.1))then
        atmnumber=42
      elseif((i.ge.97.9).and.(i.le.98.1))then
        atmnumber=43
      elseif((i.ge.109.9).and.(i.le.101.5))then
        atmnumber=44
      elseif((i.ge.102.5).and.(i.le.103.1))then
        atmnumber=45
      elseif((i.ge.105.9).and.(i.le.106.5))then
        atmnumber=46
      elseif((i.ge.107.5).and.(i.le.108.1))then
        atmnumber=47
      elseif((i.ge.111.9).and.(i.le.112.5))then
        atmnumber=48
      elseif((i.ge.114.5).and.(i.le.115.1))then
        atmnumber=49
      elseif((i.ge.118.5).and.(i.le.119.1))then
        atmnumber=50
      elseif((i.ge.121.5).and.(i.le.122.1))then
        atmnumber=51
      elseif((i.ge.127.5).and.(i.le.128.1))then
        atmnumber=52
      elseif((i.ge.126.5).and.(i.le.127.1))then
        atmnumber=53
      elseif((i.ge.130.9).and.(i.le.131.5))then
        atmnumber=54
      elseif((i.ge.132.5).and.(i.le.133.1))then
        atmnumber=55
      elseif((i.ge.136.9).and.(i.le.137.5))then
        atmnumber=56
      elseif((i.ge.138.5).and.(i.le.139.1))then
        atmnumber=57
      elseif((i.ge.139.9).and.(i.le.140.5))then
        atmnumber=58
      elseif((i.ge.140.6).and.(i.le.141.1))then
        atmnumber=59
      elseif((i.ge.144.0).and.(i.le.144.5))then
        atmnumber=60
      elseif((i.ge.144.9).and.(i.le.145.1))then
        atmnumber=61
      elseif((i.ge.150.0).and.(i.le.150.6))then
        atmnumber=62
      elseif((i.ge.151.5).and.(i.le.152.1))then
        atmnumber=63
      elseif((i.ge.156.9).and.(i.le.157.5))then
        atmnumber=64
      elseif((i.ge.158.5).and.(i.le.159.1))then
        atmnumber=65
      elseif((i.ge.162.0).and.(i.le.163.1))then
        atmnumber=66
      elseif((i.ge.164.5).and.(i.le.165.1))then
        atmnumber=67
      elseif((i.ge.166.5).and.(i.le.167.9))then
        atmnumber=68
      elseif((i.ge.168.0).and.(i.le.169.1))then
        atmnumber=69
      elseif((i.ge.172.9).and.(i.le.173.5))then
        atmnumber=70
      elseif((i.ge.174.0).and.(i.le.175.1))then
        atmnumber=71
      elseif((i.ge.178.0).and.(i.le.179.1))then
        atmnumber=72
      elseif((i.ge.180.0).and.(i.le.181.1))then
        atmnumber=73
      elseif((i.ge.183.0).and.(i.le.184.1))then
        atmnumber=74
      elseif((i.ge.185.9).and.(i.le.186.5))then
        atmnumber=75
      elseif((i.ge.189.9).and.(i.le.190.5))then
        atmnumber=76
      elseif((i.ge.191.9).and.(i.le.192.5))then
        atmnumber=77
      elseif((i.ge.194.9).and.(i.le.195.5))then
        atmnumber=78
      elseif((i.ge.196.5).and.(i.le.197.1))then
        atmnumber=79
      elseif((i.ge.200.0).and.(i.le.201.1))then
        atmnumber=80
      elseif((i.ge.203.9).and.(i.le.204.6))then
        atmnumber=81
      elseif((i.ge.206.9).and.(i.le.207.6))then
        atmnumber=82
      elseif((i.ge.208.5).and.(i.le.209.1))then
        atmnumber=83
c     Po atomic number 84 has the same mass range as Bi (83)
      elseif((i.ge.209.9).and.(i.le.210.1))then
        atmnumber=85
      elseif((i.ge.221.9).and.(i.le.222.1))then
        atmnumber=86
      elseif((i.ge.222.9).and.(i.le.223.1))then
        atmnumber=87
      elseif((i.ge.225.9).and.(i.le.226.1))then
        atmnumber=88
      elseif((i.ge.226.9).and.(i.le.227.1))then
        atmnumber=89
      elseif((i.ge.231.9).and.(i.le.232.1))then
        atmnumber=90
      elseif((i.ge.230.9).and.(i.le.231.1))then
        atmnumber=91
      elseif((i.ge.237.9).and.(i.le.238.1))then
        atmnumber=92
c     Np atomic number 93 has the same mass range as U (92)
      elseif((i.ge.243.9).and.(i.le.244.1))then
        atmnumber=94
      elseif((i.ge.242.9).and.(i.le.243.1))then
        atmnumber=95
      elseif((i.ge.246.9).and.(i.le.247.1))then
        atmnumber=96
c     Bk atomic number 97 has the same mass range as Cm (96)
      elseif((i.ge.250.9).and.(i.le.251.1))then
        atmnumber=98
      elseif((i.ge.251.9).and.(i.le.252.1))then
        atmnumber=99
      elseif((i.ge.256.9).and.(i.le.257.1))then
        atmnumber=100
      elseif((i.ge.257.9).and.(i.le.258.1))then
        atmnumber=101
      elseif((i.ge.258.9).and.(i.le.259.1))then
        atmnumber=102
      elseif((i.ge.261.9).and.(i.le.262.1))then
        atmnumber=103
      elseif((i.ge.266.9).and.(i.le.267.1))then
        atmnumber=104
      elseif((i.ge.267.9).and.(i.le.268.1))then
        atmnumber=105
      elseif((i.ge.270.9).and.(i.le.271.1))then
        atmnumber=106
      elseif((i.ge.269.9).and.(i.le.270.1))then
        atmnumber=107
      elseif((i.ge.268.9).and.(i.le.269.1))then
        atmnumber=108
      elseif((i.ge.277.9).and.(i.le.278.1))then
        atmnumber=109
      elseif((i.ge.280.9).and.(i.le.281.1))then
        atmnumber=110
c     Rg atomic number 112 has the same mass range as Ds (110)
      elseif((i.ge.284.9).and.(i.le.285.1))then
        atmnumber=112
      elseif((i.ge.285.9).and.(i.le.286.1))then
        atmnumber=113
      elseif((i.ge.288.9).and.(i.le.289.1))then
        atmnumber=114
c     Uup atomic number 115 has the same mass range as Fl (114)
      elseif((i.ge.292.9).and.(i.le.293.1))then
        atmnumber=116
      elseif((i.ge.293.9).and.(i.le.294.1))then
        atmnumber=117
c     Uuo atomic number 118 has the same mass range as Uus (117)
      endif
      return
      end function atmnumber

      subroutine timchk(ktim,time)

c***********************************************************************
c     
c     Routine for time elapsed in seconds
c     
c***********************************************************************
      implicit none

      character*12 dat,tim,zon
      integer mynode,ktim,day
      real(8) time,told,tnow
      integer info(8)

      save day

      call date_and_time(dat,tim,zon,info)
       
      if(ktim.eq.0)then

         day=info(3)
         time=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))

      else 
         told=time
         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         time=tnow-told
      endif

      return
      end subroutine timchk

      subroutine getroot(cubefile,ind)
c*********************************************************************
c     gather the filename before the .cube
c     this will enable the writing of other output files
c*********************************************************************
      implicit none
      integer i,ind
      character*100 cubefile

      do i=1,100
        if(cubefile(i:i+4).eq.'.cube')then
           ind=i-1
        endif
      enddo
      allocate(rootname(ind))
      do i=1,ind
        rootname(i)=cubefile(i:i)
      enddo
      return
      end subroutine getroot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to deal with the simulation box and atomic variables
!-----according to the PBCs    
      subroutine process_box(n_atoms,ngx,ngy,ngz)
      implicit none
      integer i_dim, j_dim, i_atom, ngx, ngy, ngz, n_atoms

      double precision a(n_dim), b(n_dim), c(n_dim), ab(n_dim),
     & r(n_dim), cross_tmp(n_dim), signo, dot, dot_1,
     & cross_tmp_1(n_dim), Box_volume

!-----Creating the real space box
      real_box_vector(1,1:n_dim) = ngx*axis_vector(1,1:n_dim)
      real_box_vector(2,1:n_dim) = ngy*axis_vector(2,1:n_dim)
      real_box_vector(3,1:n_dim) = ngz*axis_vector(3,1:n_dim)

      a(1) = real_box_vector(1,1)
      a(2) = real_box_vector(1,2)
      a(3) = real_box_vector(1,3)
      b(1) = real_box_vector(2,1)
      b(2) = real_box_vector(2,2)
      b(3) = real_box_vector(2,3)
      c(1) = real_box_vector(3,1)
      c(2) = real_box_vector(3,2)
      c(3) = real_box_vector(3,3)

      call Vect_Cross(b,c,ab)
      call Vect_Dot(a,ab,Box_volume) 
      !write(*,*) "Real box volume ", Box_volume

!-----Creating the reciprocal space box
      call Vect_Cross(a,b,ab)
      do i_dim=1, n_dim
       recip_box_vector(3,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do
      call Vect_Cross(c,a,ab)
      do i_dim=1, n_dim
       recip_box_vector(2,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do
      call Vect_Cross(b,c,ab)
      do i_dim=1, n_dim
       recip_box_vector(1,i_dim) = (2*pi/Box_volume)*ab(i_dim)
      end do

!-----Setting the origin of the real box to (0,0,0) 
      do i_atom=1, n_atoms
       do i_dim=1, n_dim
        atom_pos(i_atom,i_dim) =
     &  atom_pos(i_atom,i_dim) - axis_zero(i_dim)
        r(i_dim) = atom_pos(i_atom,i_dim)
       end do
!------Including the atoms within the box (applying pbc's)       
       call Vect_Cross(r,b,cross_tmp)
       call Vect_Dot(c,cross_tmp,dot)
       call Vect_Cross(a,b,cross_tmp_1)
       call Vect_Dot(c,cross_tmp_1,dot_1)
       signo = dot/dot_1
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) + a(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) - a(i_dim)
        end do
       end if

       call Vect_Cross(r,c,cross_tmp)
       call Vect_Dot(a,cross_tmp,dot)
       call Vect_Cross(b,c,cross_tmp_1)
       call Vect_Dot(a,cross_tmp_1,dot_1)
       signo = dot/dot_1
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) + b(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) - b(i_dim)
        end do
       end if

       call Vect_Cross(r,a,cross_tmp)
       call Vect_Dot(b,cross_tmp, dot)
       call Vect_Cross(c,a,cross_tmp_1)
       call Vect_Dot(b,cross_tmp_1,dot_1)
       signo = dot/dot_1
       if(signo.lt.0) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) + c(i_dim)
        end do
       else if(signo.gt.1) then
        do i_dim=1, n_dim
         atom_pos(i_atom,i_dim) =
     &   atom_pos(i_atom,i_dim) - c(i_dim)
        end do
       end if

       atom_pos_frac(i_atom,1) = 
     &atom_pos(i_atom,1)*recip_box_vector(1,1) + 
     &atom_pos(i_atom,2)*recip_box_vector(2,1) + 
     &atom_pos(i_atom,3)*recip_box_vector(3,1)
       atom_pos_frac(i_atom,2) = 
     &atom_pos(i_atom,1)*recip_box_vector(1,2) + 
     &atom_pos(i_atom,2)*recip_box_vector(2,2) + 
     &atom_pos(i_atom,3)*recip_box_vector(3,2)
       atom_pos_frac(i_atom,3) = 
     &atom_pos(i_atom,1)*recip_box_vector(1,3) + 
     &atom_pos(i_atom,2)*recip_box_vector(2,3) + 
     &atom_pos(i_atom,3)*recip_box_vector(3,3)
      end do

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine to compute the cross product of two vectors
      Subroutine Vect_Cross(v1,v2,v3)
      Implicit NONE

      double precision v1(3),v2(3),v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      
      End subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Function to compute the dot product of two vectors
      subroutine Vect_Dot(a,b,dot)
      Implicit NONE

      double precision a(3),b(3),dot

      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      
      end subroutine Vect_Dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----Subroutine that assigns the van der Waals radii for the
!-----elements found in the cube file according to the UFF tabulation 
      subroutine VDW_radii_array(vdw_radii,n_atoms)
      integer i,n_atoms
      double precision vdw_radii(n_atoms)
      double precision vdw_radii_file(500)

!-----Hardcore tabulation taken from the UFF with
!-----the VDW radii for the elements in the periodic table
!-----ordered according to their atomic number
!-----UNITS ARE IN BOHR
      vdw_radii_file(1) = 2.72687
      vdw_radii_file(2) = 2.23177
      vdw_radii_file(3) = 2.31586
      vdw_radii_file(4) = 2.59365
      vdw_radii_file(5) = 3.85788
      vdw_radii_file(6) = 3.63867
      vdw_radii_file(7) = 3.4582
      vdw_radii_file(8) = 3.30702
      vdw_radii_file(9) = 3.17852
      vdw_radii_file(10) = 3.06419
      vdw_radii_file(11) = 2.81853
      vdw_radii_file(12) = 2.85443
      vdw_radii_file(13) = 4.25094
      vdw_radii_file(14) = 4.05819
      vdw_radii_file(15) = 3.91835
      vdw_radii_file(16) = 3.81252
      vdw_radii_file(17) = 3.72937
      vdw_radii_file(18) = 3.65473
      vdw_radii_file(19) = 3.60182
      vdw_radii_file(20) = 3.21159
      vdw_radii_file(21) = 3.11332
      vdw_radii_file(22) = 2.99994
      vdw_radii_file(23) = 2.97065
      vdw_radii_file(24) = 2.85632
      vdw_radii_file(25) = 2.79774
      vdw_radii_file(26) = 2.75144
      vdw_radii_file(27) = 2.71365
      vdw_radii_file(28) = 2.67774
      vdw_radii_file(29) = 3.3023
      vdw_radii_file(30) = 2.61066
      vdw_radii_file(31) = 4.14133
      vdw_radii_file(32) = 4.04401
      vdw_radii_file(33) = 3.99677
      vdw_radii_file(34) = 3.97315
      vdw_radii_file(35) = 3.95803
      vdw_radii_file(36) = 3.91268
      vdw_radii_file(37) = 3.88717
      vdw_radii_file(38) = 3.44025
      vdw_radii_file(39) = 3.16057
      vdw_radii_file(40) = 2.95175
      vdw_radii_file(41) = 2.99049
      vdw_radii_file(42) = 2.88372
      vdw_radii_file(43) = 2.8327
      vdw_radii_file(44) = 2.79963
      vdw_radii_file(45) = 2.7675
      vdw_radii_file(46) = 2.73916
      vdw_radii_file(47) = 2.97443
      vdw_radii_file(48) = 2.69097
      vdw_radii_file(49) = 4.21692
      vdw_radii_file(50) = 4.14984
      vdw_radii_file(51) = 4.17629
      vdw_radii_file(52) = 4.22354
      vdw_radii_file(53) = 4.25188
      vdw_radii_file(54) = 4.16118
      vdw_radii_file(55) = 4.26795
      vdw_radii_file(56) = 3.49883
      vdw_radii_file(57) = 3.32781
      vdw_radii_file(58) = 3.35993
      vdw_radii_file(59) = 3.40718
      vdw_radii_file(60) = 3.37789
      vdw_radii_file(61) = 3.35143
      vdw_radii_file(62) = 3.32592
      vdw_radii_file(63) = 3.30041
      vdw_radii_file(64) = 3.1823
      vdw_radii_file(65) = 3.26072
      vdw_radii_file(66) = 3.23899
      vdw_radii_file(67) = 3.22104
      vdw_radii_file(68) = 3.20403
      vdw_radii_file(69) = 3.18797
      vdw_radii_file(70) = 3.17002
      vdw_radii_file(71) = 3.4393
      vdw_radii_file(72) = 2.96781
      vdw_radii_file(73) = 2.99522
      vdw_radii_file(74) = 2.89978
      vdw_radii_file(75) = 2.79113
      vdw_radii_file(76) = 2.94797
      vdw_radii_file(77) = 2.68341
      vdw_radii_file(78) = 2.60215
      vdw_radii_file(79) = 3.11143
      vdw_radii_file(80) = 2.55585
      vdw_radii_file(81) = 4.10732
      vdw_radii_file(82) = 4.06008
      vdw_radii_file(83) = 4.12905
      vdw_radii_file(84) = 4.44936
      vdw_radii_file(85) = 4.4881
      vdw_radii_file(86) = 4.50227
      vdw_radii_file(87) = 4.62983
      vdw_radii_file(88) = 3.47426
      vdw_radii_file(89) = 3.28623
      vdw_radii_file(90) = 3.20875
      vdw_radii_file(91) = 3.23521
      vdw_radii_file(92) = 3.20781
      vdw_radii_file(93) = 3.23521
      vdw_radii_file(94) = 3.23521
      vdw_radii_file(95) = 3.19458
      vdw_radii_file(96) = 3.14261
      vdw_radii_file(97) = 3.1549
      vdw_radii_file(98) = 3.13033
      vdw_radii_file(99) = 3.1171
      vdw_radii_file(100) = 3.10482
      vdw_radii_file(101) = 3.09348
      vdw_radii_file(102) = 3.06892
      vdw_radii_file(103) = 3.05758

!-----Assigning the corresponding VDW radii to the elements found in
!-----the ESP cube file
      do i=1, n_atoms
       vdw_radii(i) =  vdw_radii_file(atom_number(i))
      end do

      end subroutine

c######################################################################
c   END OF PROGRAM
c######################################################################

      end program 

