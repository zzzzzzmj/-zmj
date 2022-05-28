
!           正压原始方程模式（Barotropic primitive equation model）

!      m=20：      纬向格点数
!      n=16：      经向格点数
!      d：         网格距
!      rm：        地图放大系数
!      f：         地转参数
!      w：         工作数组
!      cla，clo：  区域中心纬度和经度
!      dt：        时间步长
!      s：         平滑系数
!      ua，ub，uc：n-1，n，n+1时间层的纬向风速
!      va，vb，vc：n-1，n，n+1时间层的经向风速
!      za，zb，zc：n-1，n，n+1时间层的位势高度
!      na：        控制12小时预报的参数
!      nb：        记录时间积分步数的参数
!      nt2=72：    判别是否积分12小时，是否该做内点平滑；
!      nt4=6：     判定是否该做边界平滑；
!      nt5=36：    判定是否该做时间平滑。 
!      zo：        为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度。
!      ni:         是否进行初始风场的静力初始化。
!                  ni=0为不进行初始化，使用读入的风场和高度场；
!                  ni=1为进行初始化，需要位势高度场做初值即可。


            program Barotropic_Model
			parameter(m=20,n=16,d=300000.0,cla=51.0,clo=118.0,dt=600.0)
            dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n),w(m,n)
			zo=2500.0
			s=0.5
            nt2=72
            nt4=6
            nt5=36
            c1=dt/2.0
            c2=dt*2.0

			write(*,*)
			write(*,*)'！！！！欢迎使用正压原始方程模式！！！！'
			write(*,*)
			write(*,*)'打开文件、读入初始场、准备输出场文件...........'

            open(11,file='D:\practice\Input\za.dat')                    ! 输入的高度场（文本）
            open(12,file='D:\practice\Input\ua.dat')                    ! 输入的u风场（文本）
            open(13,file='D:\practice\Input\va.dat')                    ! 输入的v风场（文本）

            open(17,file='D:\practice\Input\za.grd',form='binary') 		  ! 输入的高度场（二进制）	
            open(18,file='D:\practice\Input\ua.grd',form='binary') 		  ! 输入的u风场（二进制）				
            open(19,file='D:\practice\Input\va.grd',form='binary') 		  ! 输入的v风场（二进制）				
			    
            open(21,file='D:\practice\Output\rm.dat')                   ! 地图放大系数
            open(22,file='D:\practice\Output\f.dat')                    ! 地转参数

            open(23,file='D:\practice\Output\ub.dat')                   ! 静力初始化得到的u风场（文本）
			open(24,file='D:\practice\Output\ub.grd',form='binary')     ! 静力初始化得到的u风场（二进制）
            open(25,file='D:\practice\Output\vb.dat')                   ! 静力初始化得到的v风场（文本）
			open(26,file='D:\practice\Output\vb.grd',form='binary')     ! 静力初始化得到的v风场（二进制）

            open(27,file='D:\practice\Output\zc.dat')                   ! 预报的高度场（文本）
			open(28,file='D:\practice\Output\zc.grd',form='binary')     ! 预报的高度场（二进制）
			open(29,file='D:\practice\Output\uc.dat')                   ! 预报的u风场（文本）
			open(30,file='D:\practice\Output\uc.grd',form='binary')     ! 预报的u风场（二进制）
            open(31,file='D:\practice\Output\vc.dat')                   ! 预报的v风场（文本）
			open(32,file='D:\practice\Output\vc.grd',form='binary')     ! 预报的v风场（二进制）


!  读入初始资料场 
  		    read(11,'(20f6.0)')za 
  		    read(12,'(20f10.5)')ua 
  		    read(13,'(20f10.5)')va 

			write(*,*)
			write(*,*)'将初始高度场和风场写成二进制文件，便于Grads绘图......'
 			write(17)((za(i,j),i=1,m),j=1,n)   
 			write(18)((ua(i,j),i=1,m),j=1,n)   
 			write(19)((va(i,j),i=1,m),j=1,n)   

!  计算放大系数和地转参数,并写入数据文件中
			write(*,*)
			write(*,*)'计算每个格点上的地图放大系数和地转参数，并写入对应输出文件......'
            call cmf(rm,f,d,cla,m,n)
            write(21,'(20f10.5)')rm
            write(22,'(20e15.5)')f
 
			write(*,*)
			write(*,*)'静力初始化选项（0表示不进行静力初始化、1表示进行静力初始化），请输入：'
			write(*,*)
			write(*,*)'注意：如果求地转风的子程序（风场初始化）未完成，则只能输入数字0。'
			read(*,*)ni

70			write(*,*)
			write(*,*)'平滑选项：（0表示不进行平滑、1表示进行正平滑、2表示进行正逆平滑），请输入：'
			read(*,*)l
            if(l/=0 .and. l/=1 .and. l/=2 )then
                write(*,*)"平滑参数选择错误，请重新选择！"
                goto 70
            endif

  
			if(ni==1)then
!  计算地转风初值 
				write(*,*)'进行静力初始化，由风场求出高度场.....'  
				call cgw(ua,va,za,rm,f,d,m,n)
				write(23,'(20f10.5)')ua
 				write(24)((ua(i,j),i=1,m),j=1,n) 
				write(25,'(20f10.5)')va
 				write(26)((va(i,j),i=1,m),j=1,n)
			elseif(ni==0)then
				write(*,*)'不进行静力初始化使用给出的位势高度场和风场......'
			else
				write(*,*)'啊!!!输入了错误字符，请重新运行程序，输入0或1！'
				goto 90

			endif
				

!   边值传送子程序
			write(*,*)
			write(*,*)'固定边界条件赋值........'
            call tbv(ub,vb,zb,ua,va,za,m,n)
            call tbv(uc,vc,zc,ua,va,za,m,n)

!   开始预报  
			write(*,*)
			write(*,*)'开始12小时预报.......'
            do na=1,2
				nb=0
!   欧拉后差积分1小时
				do nn=1,6
					call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,dt,zo,m,n)
					call ti(ua,va,za,ub,vb,zb,ua,va,za,rm,f,d,dt,zo,m,n)
                    nb=nb+1
				enddo

!   边界平滑子程序
				call ssbp(za,w,s,m,n)
				call ssbp(ua,w,s,m,n)
				call ssbp(va,w,s,m,n)

!   前差积分半步
				call ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,c1,zo,m,n)
!   中央差积分半步
				call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
				nb=nb+1
!  数组传送子程序
				call ta(ub,vb,zb,uc,vc,zc,m,n)
!  中央差积分一步,共积分11小时
				do nn=1,66
					call ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,c2,zo,m,n)
					nb=nb+1
!  打印积分步数,na大循环步,nb小循环步
					call pv(na,nb)
!  判断是否积分12小时
					if(nb.eq.nt2) go to 80
!  判断是否做边界平滑
					if(nb/nt4*nt4.eq.nb)then
						call ssbp(zc,w,s,m,n)
						call ssbp(uc,w,s,m,n)
						call ssbp(vc,w,s,m,n)
					else	
!  判断是否做时间平滑
						if(nb.eq.nt5 .or. nb.eq.nt5+1)then
!  时间平滑子程序
                            call ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,m,n)
						else
!  数组传送,为下一轮积分做准备
							call ta(ua,va,za,ub,vb,zb,m,n)
							call ta(ub,vb,zb,uc,vc,zc,m,n)
						endif
					endif
				enddo

!  区域内点平滑
80              continue
    		    call ssip(zc,w,s,m,n,l)
				call ssip(uc,w,s,m,n,l)
				call ssip(vc,w,s,m,n,l)

!  打印积分步数
				call pv(na,nb)

!  数组传送,为后12小时的积分做准备
				call ta(ua,va,za,uc,vc,zc,m,n)
			enddo

!  存放预报结果
			write(*,*)
			write(*,*)'输出预报结果.......'
            write(27,'(20f6.0)') zc
			write(28) ((zc(i,j),i=1,m),j=1,n)    
            write(29,'(20f10.5)')uc
			write(30) ((uc(i,j),i=1,m),j=1,n) 
            write(31,'(20f10.5)')vc
            write(32) ((vc(i,j),i=1,m),j=1,n) 
			write(*,*)
			write(*,*)'！！！预报结束！！！'
90			stop
			
            end




!     computing map factors and coriolis parameter
!     rk:  圆锥常数
!     rlq: 兰勃特投影映像平面上赤道到北极点的距离
!     a:   地球半径
!     sita:标准余纬
!     psx: 区域中心余纬
!     r:   模式中心到北极的距离
            subroutine cmf(rm,f,d,cla,m,n)
            dimension rm(m,n),f(m,n)
            rk=0.7156
            rlq=11423370.0
            a=6371000.0
            conv=57.29578
            w1=2.0/rk
            sita=30.0/conv
            psx=(90.0-cla)/conv

!  计算模式中心到北极的距离r 
            cel0=a*sin(sita)/rk
            cel=(tan(psx/2.0))/(tan(sita/2.0))
            r=cel0*cel**rk

!  确定网格坐标原点在地图坐标系中的位置
            xi0=-(m-1)/2.0
            yj0=r/d+(n-1)/2.0

!  求各格点至北极点的距离rl,(xj,yi)为模式各格点在地图坐标系中的位置  
            do i=1,m
				do j=1,n
					xi=xi0+(i-1)
					yj=yj0-(j-1)
					rl=sqrt(xi**2+yj**2)*d

!  求放大系数rm和柯氏参数f
					w2=(rl/rlq)**w1
					sinl=(1.0-w2)/(1.0+w2)
					rm(i,j)=rk*rl/(a*sqrt(1.0-sinl**2))
					f(i,j)=1.4584e-4*sinl
				enddo
			enddo
            return
            end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    computing geostrophic winds
!    请同学编写地转风初值的子程序！！！应用书中（4.134）式或参见实习说明文档
            subroutine cgw(ua,va,za,rm,f,d,m,n)
            dimension ua(m,n),va(m,n),za(m,n),rm(m,n),f(m,n)
			g=9.8
			do i=1,m
				ua(i,1)=-rm(i,1)*g/f(i,1)*(za(i,1)-za(i,2))/d
				ua(i,n)=-rm(i,n)*g/f(i,n)*(za(i,n)-za(i,n-1))/d
			enddo
			do j=1,n
				va(1,j)=rm(1,j)*g/f(1,j)*(za(1,j)-za(2,j))/d
				va(m,j)=rm(m,j)*g/f(m,j)*(za(m,j)-za(m-1,j))/d
			enddo
			do i=1,m
				do j=2,n-1
					ua(i,j)=-rm(i,j)*g/f(i,j)*(za(i,j)-za(i,j))/2d
				enddo
			enddo
			do i=2,m-1
				do j=1,n
					va(i,j)=rm(1,j)*g/f(1,j)*(za(1,j)-za(2,j))/2d
				enddo
			enddo
			return
            end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     time integrations
            subroutine ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,m,n)
            dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),uc(m,n),vc(m,n),zc(m,n),rm(m,n),f(m,n)
			c=0.25/d
            m1=m-1
            n1=n-1
            do i=2,m1
				do j=2,n1
					e=-c*rm(i,j)*((ub(i+1,j)+ub(i,j))*(ub(i+1,j)-ub(i,j))+(ub(i,j)+ub(i-1,j))*(ub(i,j)-ub(i-1,j))   &
						+(vb(I,j-1)+vb(i,j))*(ub(i,j)-ub(i,j-1))+(vb(I,j)+vb(i,j+1))*(ub(i,j+1)-ub(i,j))              &
						+19.6*(zb(i+1,j)-zb(i-1,j)))+f(i,j)*vb(i,j)
					uc(i,j)=ua(i,j)+e*dt
					g=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(vb(i+1,j)-vb(i,j))+(ub(I,j)+ub(i-1,j))*(vb(i,j)-vb(i-1,j))   &
						+(vb(I,j-1)+vb(i,j))*(vb(i,j)-vb(i,j-1))+(vb(I,j)+vb(i,j+1))*(vb(i,j+1)-vb(i,j))              &
						+19.6*(zb(i,j+1)-zb(i,j-1)))-f(i,j)*ub(i,j)
					vc(i,j)=va(i,j)+g*dt
				enddo
			enddo
            do i=2,m1
				do j=2,n1
					h=-c*rm(i,j)*((ub(I+1,j)+ub(i,j))*(zb(i+1,j)-zb(i,j))+(ub(I,j)+ub(i-1,j))*(zb(i,j)-zb(i-1,j))   &
						+(vb(I,j-1)+vb(i,j))*(zb(i,j)-zb(i,j-1))+(vb(I,j)+vb(i,j+1))*(zb(i,j+1)-zb(i,j))              &
						+2.0*(zb(i,j)-zo)*(ub(i+1,j)-ub(i-1,j)+vb(i,j+1)-vb(i,j-1)))
					zc(i,j)=za(i,j)+h*dt
				enddo
			enddo
			return
            end

!     time smoothimg
            subroutine ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,m,n)
            dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n),uc(m,n),vc(m,n),zc(m,n)
            m1=m-1
            n1=n-1
            do i=2,m1
				do j=2,n1
					ub(i,j)=ub(i,j)+s*(ua(i,j)+uc(i,j)-2.0*ub(i,j))/2.0
					vb(i,j)=vb(i,j)+s*(va(i,j)+vc(i,j)-2.0*vb(i,j))/2.0
					zb(i,j)=zb(i,j)+s*(za(i,j)+zc(i,j)-2.0*zb(i,j))/2.0
				enddo
			enddo
            return
            end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     space smoothing for internal points 区域内5点平滑(正逆平滑)
!   请同学编写区域内5点平滑(正逆平滑)的子程序！！！应用书中（4.126）式或参见实习说明文档
!     l=0为不平滑，l=1为只执行正平滑，l=2为执行正逆平滑
            subroutine ssip(a,w,s,m,n,l)
            dimension a(m,n),w(m,n)     			
			s=0.5
			m1=m-1
			n1=n-1
			if(l==0)
				do i=2,m1
					do j=2,n1
						w(i,j)=a(i,j)
					enddo
				enddo

			endif
			if(l==1)
				do i=2,m1
					do j=2,n1
						w(i,j)=a(i,j)+0.25*s*(a(i+1,j)+a(i,j+1)+a(i-1,j)+a(i,j-1)-4*a(i,j))
					enddo
				enddo
				
			endif
			if(l==2)
				do 
			
  		    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!     space smoothing for boundary points 边界九点平滑
            subroutine ssbp(a,w,s,m,n)
            dimension a(m,n),w(m,n)
            m1=m-1
            m3=m-3
            n1=n-1
            n2=n-2
            n3=n-3
            do i=2,m1
				do j=2,n1,n3
					w(i,j)=a(i,j)+0.5*s*(1.0-s)*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))    &
						    +0.25*s*s*(a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
				enddo
			enddo
            do i=2,m1,m3
				do j=3,n2
					w(i,j)=a(i,j)+0.5*s*(1.0-s)*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.0*a(i,j))    &
							+0.25*s*s*(a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1)-4.0*a(i,j))
				enddo
			enddo
            do i=2,m1
				do j=2,n1,n3
					a(i,j)=w(i,j)
				enddo
			enddo
            do i=2,m1,m3
				do j=3,n2
					a(i,j)=w(i,j)
				enddo
			enddo
			return
            end


!     transmiting arrays  数组传送
            subroutine ta(ua,va,za,ub,vb,zb,m,n)
            dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n)
            do i=1,m
				do j=1,n
					ua(i,j)=ub(i,j)
					va(i,j)=vb(i,j)
					za(i,j)=zb(i,j)
				enddo
			enddo
            return
            end

!     transmiting boundary valaus  赋固定边界值
            subroutine tbv(ua,va,za,ub,vb,zb,m,n)
            dimension ua(m,n),va(m,n),za(m,n),ub(m,n),vb(m,n),zb(m,n)
            m1=m-1
            n1=n-1
            do i=1,m
				do j=1,n,n1
					ua(i,j)=ub(i,j)
					va(i,j)=vb(i,j)
					za(i,j)=zb(i,j)
				enddo
			enddo
            do i=1,m,m1
				do j=1,n
					ua(i,j)=ub(i,j)
					va(i,j)=vb(i,j)
					za(i,j)=zb(i,j)
				enddo
			enddo
            return
            end

!     printing variables  打印积分步数
            subroutine pv(na,nb)
            write(*,'(5x,3hna=,i3,5x,3hnb=,i2/)')na,nb
            return
            end 