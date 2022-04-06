      module neib
!     ! for given plat,pdiff rmax, return dislist
      integer,protected:: ng
      real(8),allocatable,protected:: dislist(:)
      contains
      subroutine neighbourpair(plat,rbas21,rmax)
!! For given plat and rbas21=ratom2-ratom1 (cartesian),
!! we give distance lists, within rmax.
      
!! -----------------------------------------------------------------------------------      
!! Show possible distance of lattice as  |rbas21- (n1*plat1 + n2*plat2 + n3*plat3)|< rmax
      implicit none
      real(8):: plat(3,3),bohr,xv,alat,pi=4*atan(1d0),qlat(3,3),rmax
     &,pwgmax,qdiff(3)
     &,dum,platx(3,3),pwgmin,gmax0,rvec(3),distance,rbas21(3),rmax2,dis2
      integer,allocatable:: igv2x(:,:),kv_iv(:)
      integer::ig,napw,nn1,nn2,nn3,mshlst,nmin(3),nmax(3),n1,n2,n3,ndatx
c      print *, plat(1:3,1),plat(1:3,2),plat(1:3,3),rbas21(1:3),rmax
c     PLAT(1:3,1)=[0.5, 0.5, 1.0]
c      PLAT(1:3,2)=[0.5, 1.0, 0.5]
c      PLAT(1:3,3)=[1.0, 0.5, 0.5]
c      rbas1=[0,0,0]       !cartesian atom1
c      rbas2=[.5,.5,0.0]    !cartesian atom2
!
      rmax2=rmax**2
      call gvlstnx(plat(1,1),plat(1,2),plat(1,3),rbas21,rmax,
     &     nmin(1),nmax(1))
      call gvlstnx(plat(1,2),plat(1,3),plat(1,1),rbas21,rmax,
     &     nmin(2),nmax(2))
      call gvlstnx(plat(1,3),plat(1,1),plat(1,2),rbas21,rmax,
     &     nmin(3),nmax(3))
!! For given rbas21(1:3) (cartesian),
c$$$      print *,'rbas21 =',rbas21
c$$$      print *,'plat1 =',plat(:,1)
c$$$      print *,'plat2 =',plat(:,2)
c$$$      print *,'plat3 =',plat(:,3)
c      print *,'nmin1 mnax1 =',nmin(1),nmax(1),'nmin2 mnax2 =',nmin(2),nmax(2),'nmin3 mnax3 =',nmin(3),nmax(3)
      ndatx=(nmax(1)-nmin(1)+1)*(nmax(2)-nmin(2)+1)*(nmax(3)-nmin(3)+1)
      if(allocated(dislist)) deallocate(dislist)
      allocate(dislist(ndatx))
      ig=0
      do n1 = nmin(1),nmax(1)
      do n2 = nmin(2),nmax(2)
      do n3 = nmin(3),nmax(3)
         rvec = rbas21 - (plat(:,1)*n1 + plat(:,2)*n2 + plat(:,3)*n3)
         dis2=sum(rvec*rvec)
         if( dis2<= rmax2) then
c            write(6,"(f10.5)") sqrt(dis2)
            ig=ig+1
            dislist(ig)=sqrt(dis2)
         endif
      enddo
      enddo
      enddo
      ng=ig
      end subroutine neighbourpair
      end module neib
!!------------------------------------------------------------

      
!     !
      subroutine gvlstnx(q0,q1,q2,qp,rmax, nmin,nmax)
!! Multiples of primitive vectors.
!! we look for allow n; n is nmin<n< nmax, satisfying  |qp+n1*q1+n2*q2+n3|< rmax
C     ----------------------------------------------------------------------
Ci   q0    : first  primitive vector
Ci   q1,q2 : other two primitive vectors.
Ci   qp    :  added to G vectors: sphere is centered G=qp.
Ci  rmax  : cutoff Radius
Co Outputs
Co   nmin  :search for lattice vectors limited to (nmin..nmax)*q0
Co   nmax  :search for lattice vectors limited to (nmin..nmax)*q0
!! NOTE: n1=constant makes a plane spanned by q1 and q2. We calculate distance to the plane 
!! (n1,n2,n3) on the plane has a distane larger than the distance.
      implicit none
      integer nmin,nmax,nt1,nt2
      real(8):: qperp(3),ddot,qqperp,q0(3),q1(3),q2(3),qp(3)
     &     ,rmax,distance1,distanceqp
C ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
      qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2) !external product qprep=q1 x q2
      qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3) !this is normal to the plane for given n1
      qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
      qperp = qperp/sqrt(sum(qperp*qperp))  !normalize vector. sum(a,b) is innner product
      distanceqp = sum(qp*qperp) !signed distance to qp from origin
      distance1=   sum(q0*qperp) !signed distance to plane for n1=1 from origin
!! we search n1 for |distance1*n1 - distanceqp|<rmax
!! fortran int gives an integer closer to zero for given float number.
      Qqperp = sum(q0*qperp)  !
      nmax = int(( rmax+distanceqp)/distance1)
      nmin = int((-rmax+distanceqp)/distance1)
      if(nmax<nmin) then
         nt1=nmin
         nt2=nmax
         nmax=nt1
         nmin=nt2
      endif
      end
      
!!----------------------------------------------------------      
      program main
      use neib,only: ng, dislist, Neighbourpair
      implicit none
      integer::u,v,i,npair,ix,id1,id2
      real(8):: plat(3,3),rbas12(3),rmax
      character*120:: aaa,bbb,filename
      read(*,*) filename
      open(newunit=u, file=trim(filename)//'_inp.dat', action='read')
      read(u,*) plat(1:3,1),plat(1:3,2),plat(1:3,3)
c      write(6,*) plat
      read(u,*) rmax  
      open(newunit=v, file=trim(filename)//'_out.dat', action='write')
      ix=0
      do
        ix=ix+1
        read(u,*,end=1019) id1,id2,rbas12(:)
        call neighbourpair(plat,rbas12,rmax)
        write(v,"(2i10,1000f10.5)") id1,id2,dislist(1:ng)
      enddo   
 1019 close(u)
      close(v)
      end

      

