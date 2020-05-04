c   AUTHOR: Manuel Almagro (malmriv@correo.ugr.es)
c   DESCRIPTION. This is a molecular dynamics engine that uses the velocity
c   Verlet algorithm to simulate a galaxy. (The first hundred-or-so lines are
c   a random number generator provided by professors at the University of Gra-
c   nada). The actual program starts after this module, and then you'll find
c   some functions and subroutines used for ad-hoc array manipulation. The
c   technical details will be told in an upcoming assignment that will be uplo-
c   aded to this GitHub repository.

        module randomnumber
        integer*8::ip,iq,is,np,nbit,ic
        parameter(ip=1279)
        parameter(iq=418)
        parameter(is=ip-iq)
        parameter(np=14)
        parameter(nbit=31)
        integer*8,dimension(1:ip)::ix
        save ic,ix
        contains

        real*8 function gaus()
        real*8::v1,v2,fac,r,gset
        integer*8::iset
        save iset,gset

        if(iset.eq.0)then
          r=2.d0
          do while(r>1.or.r.eq.0.d0)
            v1=2.d0*dran_u()-1.d0
            v2=2.d0*dran_u()-1.d0
            r=v1*v1+v2*v2
          enddo
          fac=dsqrt(-2.d0*dlog(r)/r)
          gset=v1*fac
          gaus=v2*fac
          iset=1
        else
          gaus=gset
          iset=0
        endif
        end function gaus

        subroutine dran_ini(iseed0)
          integer*8::m,np1,nn,nn1,i,j
          real*8::dseed,p,t,x
          dseed=iseed0
          do i=1,ip
            ix(i)=0
            do j=0,nbit-1
              if(rand_xx(dseed).lt.0.5d0) ix(i)=ibset(ix(i),j)
            enddo
          enddo
          ic=0
        end subroutine dran_ini

        subroutine dran_read(iunit)
          integer*8::i
          read(iunit,*)ic
          read(iunit,*)(ix(i),i=1,ip)
        end subroutine dran_read

        subroutine dran_write(iunit)
          integer*8::i
          write(iunit,*) ic
          write(iunit,*) (ix(i),i=1,ip)
        end subroutine dran_write

        integer*8 function i_dran(n)
        integer*8::i_ran,n
        ic=ic+1
        if(ic.gt.ip) ic=1
        if(ic.gt.iq)then
          ix(ic)=ieor(ix(ic),ix(ic-iq))
        else
          ix(ic)=ieor(ix(ic),ix(ic+is))
        endif
        i_ran=ix(ic)
        if(n.gt.0)i_dran=mod(i_ran,n)+1
        end function i_dran

        real*8 function dran_u()
        real*8::rmax
        parameter (rmax=2147483647.0)
        ic=ic+1
        if(ic.gt.ip) ic=1
        if(ic.gt.iq)then
          ix(ic)=ieor(ix(ic),ix(ic-iq))
        else
          ix(ic)=ieor(ix(ic),ix(ic+is))
        endif
        dran_u=dble(ix(ic))/rmax
        end function dran_u

        real*8 function rand_xx(dseed)
        real*8:: a,c,xm,rm,dseed
        parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
        dseed=mod(dseed*a+c,xm)
        rand_xx=dseed*rm
        end function rand_xx

        end module randomnumber

c   IMPORTANT. The actual program starts here.

        program galaxy
          use randomnumber
          implicit none
          real*8, dimension (:,:), allocatable :: x,y,ax,ay,KE,UE
          real*8, dimension (:,:), allocatable :: vx,vxcol,vy,vycol
          real*16, dimension (:), allocatable :: m, R
          integer, dimension(:,:), allocatable :: vlist
          real*8 h, t, tmax, pos1, pos2, dp, dist, vradius
          real*8 calc_x, calc_v, calc_a, coll, total_KE, total_UE
          integer i, j, k, l, n_bodies, N, nonzero, randseed(1:8)
          real*8 start, finish, inertia
          integer*8 seed1byte

          call cpu_time(start)
          n_bodies = 1517
          N = 200
          tmax = 10.d0 !1 unit ~ current age of the sun
          h = tmax/float(N)
          vradius = 0.4 !In normalised distance units

c   Allocate memory for each vector
          allocate(x(1:n_bodies,1:N),y(1:n_bodies,1:N))
          allocate(vx(1:n_bodies,1:N),vy(1:n_bodies,1:N))
          allocate(vxcol(1:n_bodies,1:N),vycol(1:n_bodies,1:N))
          allocate(ax(1:n_bodies,1:N),ay(1:n_bodies,1:N))
          allocate(KE(1:n_bodies,1:N),UE(1:n_bodies,1:N))
          allocate(m(1:n_bodies),R(1:n_bodies))
          allocate(vlist(1:n_bodies,1:n_bodies))

c   Initialize as zero the vectors whose definition implies addition
          UE = 0.d0
          KE = 0.d0
          t = 0.d0
          ax = 0.d0
          ay = 0.d0

c   Data is generated and normalised by an R script.
          open(100,file="./dataset.txt",status='old')
          do i=1,n_bodies
            read(100,*) m(i), x(i,1), y(i,1), vx(i,1), vy(i,1), R(i)
          end do
          close(100)

c   We open different files for different magnitudes. This may seem
c   inefficient for a small number of bodies, but a simulation with
c   hundreds of bodies, each one evaluated thousands of times,
c   would quickly yield a very unmanageable .txt file if we were to
c   store everything in one go.
          call system("mkdir -p ./results")
          open(10,file="./results/x.txt")
          open(20,file="./results/y.txt")
          open(30,file="./results/vx.txt")
          open(40,file="./results/vy.txt")
          open(50,file="./results/ax.txt")
          open(60,file="./results/ay.txt")
          open(70,file="./results/t.txt")
          open(80,file="./results/energy.txt")
          open(110,file="./results/blackholemass.txt")
          open(120,file="./results/inertia.txt")

c   Evaluate accelerations to set initial conditions
          i = 1
          do j=1,n_bodies !j = body whose motion is being computed
            do k=1,n_bodies !k = its interactions with the k-th body
            if(j .ne. k) then !Avoid self-interactions
            ax(j,i) = ax(j,i) + calc_a(m(k),x(j,i),x(k,i),y(j,i),y(k,i))
            ay(j,i) = ay(j,i) + calc_a(m(k),y(j,i),y(k,i),x(j,i),x(k,i))
            end if
            end do
          end do

c   Compute initial Verlet list (only every few steps)
          vlist = 0
          do j=1,n_bodies
            do k=1,n_bodies
              if(j.eq.k) cycle
              dist = (x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0
              if(dist .lt. vradius**2.0) then
                vlist(j,k) = k
              end if
            end do
            call quicksort(vlist(j,1:n_bodies),1,n_bodies)
            call reverse(vlist(j,1:n_bodies),n_bodies)
            vlist(j,1) = 1
          end do

c   Save the initial time
          write(70,*) t

c   Initiate random number generator using current time as seed
          call date_and_time(values=randseed)
          call dran_ini(abs(product(randseed)))

c   Iterative process computing motion magnitudes for each step
          do i=2,N !Step number. Indexes mean the same as before.

            do j=1,n_bodies
c   Compute positions
              x(j,i) = calc_x(h,x(j,i-1),vx(j,i-1),ax(j,i-1))
              y(j,i) = calc_x(h,y(j,i-1),vy(j,i-1),ay(j,i-1))
            end do

c   Compute Verlet list ever 4 steps, best to determine later
            if(modulo(i,4) .eq. 0) then
              vlist = 0
              do j=1,n_bodies
                do k=1,n_bodies
                  if(j.eq.k) cycle
                  dist = (x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0
                  if(dist .lt. vradius**2.0) then
                    vlist(j,k) = k
                  end if
                end do
                call quicksort(vlist(j,1:n_bodies),1,n_bodies)
                call reverse(vlist(j,1:n_bodies),n_bodies)
              vlist(j,1) = 1
              end do
            end if

c   Compute accelerations (except if there are no neighbours: stray bodies)
            do j=1,n_bodies
              if(nonzero(vlist(j,1:n_bodies),n_bodies) .eq. 0) cycle
              do k=1,nonzero(vlist(j,1:n_bodies),n_bodies)
                l = vlist(j,k)
                if(l .eq. j) cycle !Avoid possible self-interaction
                ax(j,i) = ax(j,i) + calc_a(m(l),x(j,i),x(l,i),
     &          y(j,i),y(l,i))
                ay(j,i) = ay(j,i) + calc_a(m(l),y(j,i),y(l,i),
     &          x(j,i),x(l,i))
              end do
            end do

c   Compute velocities
            do j=1,n_bodies
              vx(j,i) = calc_v(h,vx(j,i-1),ax(j,i-1),ax(j,i))
              vy(j,i) = calc_v(h,vy(j,i-1),ay(j,i-1),ay(j,i))
            end do

c   Save velocities to compute collisions
            vxcol = vx
            vycol = vy

c   Compute velocities in case a collision has taken place.
c   Define the necessary parameters to check for collisions
            do j=1,(n_bodies-1)
              do k=(j+1),n_bodies
                dist = dsqrt((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
                coll = dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),x(k,i)-x(j,i),
     &          y(k,i)-y(j,i))

c   The central black hole swallows anything that comes too close and
c   a new body with fitting magnitudes is generated.
                if(dist .le. (R(k)+R(j)) .and. j.eq.1 .or. k.eq.1) then
                  l = k
                  if(m(j) .lt. m(k)) l = j
                  m(1) = m(1) + m(l)
                  x(l,i) = (-1.0)**i_dran(seed1byte)*dran_u()*sqrt(4.5) !Galactic radius
                  y(l,i) = (-1.0)**i_dran(seed1byte)*dran_u()*sqrt(4.5)
                  vx(l,i) = y(l,i)/(x(l,i)**2.d0+y(l,i)**2.d0)
                  vy(l,i) = -x(l,i)/(x(l,i)**2.d0+y(l,i)**2.d0)
                  m(l) = 2.621810E-7+(-1)**i_dran(seed1byte)*
     &            2.621810E-7*dran_u()
                  write(110,*) t,m(1)
                end if

                if((dist .le. (R(k)+R(j))) .and. (coll .le. 0.d0)) then
c   This horrible mess computes after-collision velocities. It looks
c   awful but I promise it does what it needs to do and it works.

c   First body's velocity
                vxcol(j,i) = vx(j,i) -2.*m(k)/(m(j)+m(k))
     &          * dp(vx(j,i)-vx(k,i),vy(j,i)-vy(k,i),
     &          x(j,i)-x(k,i),y(j,i)-y(k,i))
     &          / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &          * (x(j,i)-x(k,i))

                vycol(j,i) = vy(j,i) -2.*m(k)/(m(j)+m(k))
     &          * dp(vx(j,i)-vx(k,i),vy(j,i)-vy(k,i),
     &          x(j,i)-x(k,i),y(j,i)-y(k,i))
     &          / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &          * (y(j,i)-y(k,i))

c   Second body's velocity
                vxcol(k,i) = vx(k,i) -2.*m(j)/(m(j)+m(k))
     &          * (dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),
     &          x(k,i)-x(j,i),y(k,i)-y(j,i))
     &          / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0))
     &          * (x(k,i)-x(j,i))

                vycol(k,i) = vy(k,i) -2.*m(j)/(m(j)+m(k))
     &          *  dp(vx(k,i)-vx(j,i),vy(k,i)-vy(j,i),
     &          x(k,i)-x(j,i),y(k,i)-y(j,i))
     &          / ((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
     &          * (y(k,i)-y(j,i))
                end if
              end do
            end do

c   Set new velocities
            vx = vxcol
            vy = vycol

c   Compute kinetic and potential energies
            do j = 1,n_bodies
              KE(j,i) = 0.5d0*m(j)*(vx(j,i)**2.d0+vy(j,i)**2.d0)
            end do

            do j=1,n_bodies
              do k=j+1,n_bodies
                dist = dsqrt((x(j,i)-x(k,i))**2.0+(y(j,i)-y(k,i))**2.0)
                UE(j,i) = UE(j,i) - (m(j)*m(k))/dist
              end do
            end do

c   Compute inertia for current state
            inertia = 0
            do j = 1,n_bodies
              inertia = inertia+(x(j,i)**2.d0+y(j,i)**2.d0)*m(j)
            end do
            write(120,*) inertia


c   Update time, save time, and save black hole mass
            t=t+h
            write(70,*) t
            if(modulo(i,N).eq.10) write(*,*) "Completed %:",i/N
          end do

c   Save kinematic computations. Remember:
c   i-th row = i-th iteration, j-th column = j-th object
          do i=1,N
            write(10,*) (x(j,i),j=1,n_bodies)
            write(20,*) (y(j,i),j=1,n_bodies)
            write(30,*) (vx(j,i),j=1,n_bodies)
            write(40,*) (vy(j,i),j=1,n_bodies)
            write(50,*) (ax(j,i),j=1,n_bodies)
            write(60,*) (ay(j,i),j=1,n_bodies)
          end do

c   Save energy computations.
          do i=2,N
            total_KE = 0.d0
            total_UE = 0.d0
            do j=1,n_bodies
              total_KE = total_KE + KE(j,i)
              total_UE = total_UE + UE(j,i)
            end do
            write(80,*) i,total_KE,total_UE
          end do

c   Close every file
          do i=10,120,10
            close(i)
          end do

c   Measure runtime and finish program
          call cpu_time(finish)
          write(*,*) "Total runtime (in seconds):",finish-start
          write(*,*) "Data saved in ./results/"
        end program galaxy

c   These functions evaluate angular position, velocity and accelera-
c   tion. Defining them here will make this program more adaptable.

        function calc_x(h,x,v,a) result(ans)
          implicit none
          real*8 h, x, v, a, ans
          ans = x + h*v + ((h**2.0)/2.0)*a
          return
        end function calc_x

        function calc_v(h,v,a1,a2) result(ans)
          implicit none
          real*8 h, v, a1, a2, ans
          ans = v + (h/2.0)*(a1+a2)
          return
        end function calc_v

        function calc_a(m,x1,x2,y1,y2) result(ans)
          implicit none
          real*8 x1,x2,y1,y2,dist,ans
          real*16 m
          dist = dsqrt((x1-x2)**2.0+(y1-y2)**2.0)**3.0
          ans =  -m*(x1-x2)/dist
          return
        end function calc_a

c     Dot product function
        function dp(u1,u2,v1,v2) result(ans)
          implicit none
          real*8 u1,u2,v1,v2,ans
          ans = u1*v1+u2*v2
          return
        end function dp

c     This is a quicksort algorithm.
        recursive subroutine quicksort(a, first, last)
          implicit none
          integer a(*), first, last, i, j
          real*8 x, t

          x = a( (first+last) / 2 )
          i = first
          j = last
          do
            do while (a(i) .lt. x)
              i=i+1
            end do
            do while (x .lt. a(j))
              j=j-1
            end do
            if (i .ge. j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            i=i+1
            j=j-1
          end do
          if (first < i-1) call quicksort(a, first, i-1)
          if (j+1 < last)  call quicksort(a, j+1, last)
        end subroutine quicksort

c     This function reverses the order of an integer array. It's used in
c     to put the non-zero elements of the Verlet list first.
        subroutine reverse(a,length)
          implicit none
          integer a(*), i, length, aux
          do i=1,floor(length/2.0)
            aux = a(i)
            a(i) = a(length-i+1)
            a(length-i+1) = aux
          end do
        end subroutine reverse

c     This function returns the number of non-zero elements in an integer
c     array. It's used to retrieve relevant information from the Verlet list.
        function nonzero(a,length) result(count)
          implicit none
          integer a(*), length, i, count
          count = 0
          do i=1,length
            if(a(i) .ne. 0) count = count+1
          end do
        end function nonzero
