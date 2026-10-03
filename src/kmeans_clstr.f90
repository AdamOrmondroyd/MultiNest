! soft k-means clustering with k=2
! D.J.C. MacKay, Information Theory, Inference & Learning Algorithms, 2003, p.304
! Aug 2006

module kmeans_clstr
  use RandomNS
  use utils1
  implicit none
      

contains
  
!----------------------------------------------------------------------
  
  !simple k-means
  subroutine kmeans3(k,pt,npt,numdim,means,cluster,min_pt)
	implicit none
 
    	integer k,numdim,npt,i1,i2,i3
    	double precision pt(numdim,npt)
    	double precision means(k,numdim),dis(min_pt,2)
    	integer cluster(npt),old_cluster(npt),r(npt,k),totR(k)
    	double precision temp,dist
    	integer i,j,x,nochg,scrap(k),min_pt
    	logical clstrd,flag
    	double precision urv,d1

    	if(k>npt/min_pt+1) k=npt/min_pt+1
	
	if(k==1) then
		do i=1,numdim
			means(1,i)=sum(pt(i,1:npt))/dble(npt)
		enddo
		cluster=1
		return
	endif
	
    	clstrd=.false.
    	nochg=0
    
!    call ranmarns(urv)
!    x=int(dble(npt)*urv)+1
!    means(1,:)=pt(:,x)
!    temp=-1.
!    do i=2,k
!    	lmean(:)=sum(means(1:i-1,:))/dble(i-1)
!    	do j=1,npt
!      	dist=sum((lmean(:)-pt(:,j))**2.)
!            if(dist>temp) then
!			temp=dist
!                  x=j
!            end if
!	end do
!	means(i,:)=pt(:,x)
!    end do

!    lmean(:)=means(1,:)
!    do i=1,k/2+1
!    	do
!      	flag=.false.
!    		call ranmarns(urv)
!    		x=int(dble(npt)*urv)+1
!      	do j=1,i-1
!      		if(x==scrap(j)) then
!            		flag=.true.
!                  	exit
!			end if
!		end do
!            if(.not.flag) exit
!	end do
!      scrap(i)=x
!    	means(i*2-1,1:numdim)=pt(1:numdim,x)
!      if(i*2>k) exit
!    	means(i*2,1:numdim)=2.*lmean(1:numdim)-pt(1:numdim,x)
!    end do

    	!choose random points as starting positions
    	do i=1,k
    		do
      			flag=.false.
            		urv=ranmarns(0)
    			x=int(dble(npt)*urv)+1
      			do j=1,i-1
      				if(x==scrap(j)) then
                  			flag=.true.
                        		exit
				endif
			enddo
            		if(flag) then
            			cycle
            		else
            			scrap(i)=x
                  		exit
			endif
		enddo
      		means(i,:)=pt(:,x)
    	enddo
    
    
    	old_cluster=0
    
    	do
    		do i=1,npt
    			temp=huge(1.d0)
			do j=1,k
				dist=sum((means(j,:)-pt(:,i))**2.)
                  		if(dist<temp) then
                  			temp=dist
                        		x=j
                  		endif
			enddo
            		r(i,:)=0
            		r(i,x)=1
			cluster(i)=x
		enddo
	
		clstrd=.true.
		do i=1,npt
			if(old_cluster(i)/=cluster(i)) then
				clstrd=.false.
				exit
			endif
		enddo
      
		if(clstrd) then
			!check if all the clusters have more than min_pt points
			do i=1,k
				if(totR(i)<min_pt) then
					dis=1.d99
					i1=min_pt-totR(i)
					do j=1,npt
						if(cluster(j)/=i .and. totR(cluster(j))>min_pt) then
							d1=sum((means(i,:)-pt(:,j))**2.)
							i3=0
							do i2=i1,1,-1
								if(d1<dis(i2,1)) then
									i3=i2
								else	
									exit
								endif
							enddo
							if(i3/=0) then
								dis(i3+1:i1,:)=dis(i3:i1-1,:)
								dis(i3,1)=d1
								dis(i3,2)=dble(j)
							endif
						endif
					enddo
					do j=1,i1
						i3=int(dis(j,2))
						i2=cluster(i3)
						cluster(i3)=i
						totR(i)=totR(i)+1
						totR(i2)=totR(i2)-1
						r(i3,:)=0
            					r(i3,i)=1
					enddo
				endif
			enddo
			do i=1,k
				!update means
				do j=1,numdim
					means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
				enddo
			enddo
			exit
		endif
      
		old_cluster=cluster
	
		do i=1,k
			totR(i)=sum(r(:,i))
            		if(totR(i)==0) cycle
		
			!update means
			do j=1,numdim
				means(i,j)=sum(r(:,i)*pt(j,:))/totR(i)
			enddo
		enddo
    	enddo
  
	
  end subroutine kmeans3
 
!---------------------------------------------------------------------- 

  !Dinosaur clustering
  function Dmeans(k,pt,npt,like,ndim,cluster,min_pt,nptk,meank,covmatk,invcovk,tmatk,evalk,eveck,kfack,effk, &
  detcovk,volk,pVol,fVol,cSwitch,nCdim)
    implicit none
    
    !input variables
    integer k !no. of clusters required
    integer npt !no. of points
    double precision like(npt) !scaled log-like
    integer ndim !dimensionality
    double precision pt(ndim,npt) !points
    integer min_pt !min no. of points allowed in a cluster
    double precision pVol !prior volume
    logical cSwitch
    integer nCdim
    
    !input/output variables
    integer nptk(k) !input: no. of points in each of k-1 clusters, output: no. of points in each of k clusters
    integer cluster(npt) !input: cluster membership of each point for k-1 clusters, 
    		    !output: cluster membership of each point for k clusters
    double precision meank(k,ndim) !input: means of k-1 clusters, output: means of k clusters
    double precision covmatk(k,ndim,ndim) !input: covmat of k-1 clusters, output: covmat of k clusters
    double precision invcovk(k,ndim,ndim) !input: invcov of k-1 clusters, output: invcov of k clusters
    double precision tmatk(k,ndim,ndim) !input: tmat of k-1 clusters, output: tmat of k clusters
    double precision evalk(k,ndim) !input: eval of k-1 clusters, output: eval of k clusters
    double precision eveck(k,ndim,ndim) !input: evec of k-1 clusters, output: evec of k clusters
    double precision kfack(k) !input: kfac of k-1 clusters, output: kfac of k clusters
    double precision effk(k) !input: eff of k-1 clusters, output: eff of k clusters
    double precision detcovk(k) !input: detcov of k-1 clusters, output: detcov of k clusters
    double precision volk(k) !input: vol of k-1 clusters, output: vol of k clusters
    double precision fVol !volume of the father ellipsoid
    
    !output variables
    integer Dmeans !0 if found a better partition, 1 if found a partition but not better, 2 if couldn't partition
    
    !work variables
    integer i,j,x,i1,i2,indx(1),count,gcount,cls(2)
    double precision, allocatable :: h(:,:),mdis(:,:),ptk(:,:),mu_tmp(:,:)
    double precision d1,d2
    logical flag
    logical, allocatable :: doCal(:), broken(:)
    !ellipsoid properties
    integer, allocatable :: nptx(:), clusterx(:),cluster2(:)
    double precision, allocatable :: meanx(:,:),covmatx(:,:,:),invcovx(:,:,:),tmatx(:,:,:),evalx(:,:)
    double precision, allocatable :: evecx(:,:,:),kfacx(:),effx(:),detcovx(:),volx(:)
    
    
    !sanity check
    if(npt<min_pt*k) then
    	Dmeans=2
	return
    endif
    
    
    allocate( h(k,npt), mdis(k,npt), ptk(ndim+1,npt), mu_tmp(2,ndim+1), meanx(k,ndim), covmatx(k,ndim,ndim), &
    invcovx(k,ndim,ndim), tmatx(k,ndim,ndim), evalx(k,ndim), evecx(k,ndim,ndim), kfacx(k),effx(k), detcovx(k), volx(k) )
    allocate( doCal(k), broken(k) )
    allocate( nptx(k), clusterx(npt), cluster2(npt) )
    
    !initialization
    nptx=nptk
    clusterx=cluster
    meanx=meank
    covmatx=covmatk
    invcovx=invcovk
    tmatx=tmatk
    evalx=evalk
    evecx=eveck
    kfacx=kfack
    effx=effk
    detcovx=detcovk
    volx=volk
    gcount=1
    broken=.false.
    
    do
    	flag=.false.
    	!first find the ellipsoid with most fractional wastage
    	d2=-1.d99
    	do i=1,k-1
    		if(broken(i)) cycle
    		d1=(volx(i)-pVol*nptx(i)/npt)/(pVol*nptx(i)/npt)
		if(d1>d2) then
			flag=.true.
			d2=d1
			j=i
		endif
    	enddo
    
    	!if none of the ellipsoids can be split then return
    	if(.not.flag) exit
	
    	broken(j)=.true.
	
	!split ellipsoid j
	i1=0
	do i=1,npt
		if(clusterx(i)==j) then
			i1=i1+1
			if(cSwitch) then
				ptk(1:nCdim,i1)=pt(1:nCdim,i)
				ptk(nCdim+1,i1)=like(i)
			else
				ptk(1:ndim,i1)=pt(1:ndim,i)
				ptk(ndim+1,i1)=like(i)
			endif	
		endif
	enddo
	
	i2=2
!	if(cSwitch) then
!		call Kmeans3(i2,ptk(1:nCdim+1,1:i1),i1,nCdim+1,mu_tmp(:,1:nCdim+1),cluster2(1:i1),min_pt)
!	else
!		call Kmeans3(i2,ptk(1:ndim+1,1:i1),i1,ndim+1,mu_tmp(:,1:ndim+1),cluster2(1:i1),min_pt)
!	endif
	if(cSwitch) then
		call Kmeans3(i2,ptk(1:nCdim,1:i1),i1,nCdim,mu_tmp(:,1:nCdim),cluster2(1:i1),min_pt)
	else
		call Kmeans3(i2,ptk(1:ndim,1:i1),i1,ndim,mu_tmp(:,1:ndim),cluster2(1:i1),min_pt)
	endif
	
	!separate out the points
	nptx(j)=0
	nptx(k)=0
	i2=1
	do i=1,i1
		do x=i2,npt
			if(clusterx(x)==j) then
				if(cluster2(i)==1) then
					clusterx(x)=j
					nptx(j)=nptx(j)+1
				else
					clusterx(x)=k
					nptx(k)=nptx(k)+1
				endif
				i2=x+1
				exit
			endif
		enddo
	enddo
	
	cls(1)=j
	cls(2)=k
	do i2=1,2
		x=0
		!separate out the points
		do i=1,npt
			if(clusterx(i)==cls(i2)) then
				x=x+1
				ptk(1:ndim,x)=pt(1:ndim,i)
			endif
		enddo
		!min volume this ellipsoid should occupy
		d1=pVol*x/npt
		!calculate the ellipsoid properties
		call CalcEllProp(x,ndim,ptk(1:ndim,1:x),meanx(cls(i2),:),covmatx(cls(i2),:,:), &
		invcovx(cls(i2),:,:),tmatx(cls(i2),:,:),evecx(cls(i2),:,:),evalx(cls(i2),:),detcovx(cls(i2)),kfacx(cls(i2)), &
		effx(cls(i2)),volx(cls(i2)),d1,.false.)
	enddo
	
	!find if it's a better partition
!	if(nptx(j)<min_pt .or. nptx(k)<min_pt) then
!		count=0
!		Dmeans=2
!	else
		count=1
		if(gcount==1) then
			!save the best partition
			nptk=nptx
    			cluster=clusterx
    			meank=meanx
    			covmatk=covmatx
   			invcovk=invcovx
   			tmatk=tmatx
    			evalk=evalx
    			eveck=evecx
    			kfack=kfacx
   			effk=effx
    			detcovk=detcovx
    			volk=volx
		endif
		if(sum(volx(1:k))<fVol) then
			Dmeans=0
		else
			Dmeans=1
		endif
!	endif
    
    	doCal=.true.
		
	!calculate the distance measure
!    	do i=1,k
!		do j=1,npt
!			mdis(i,j)=MahaDis(ndim,pt(:,j),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
!		enddo
!    	enddo
    
    	do
    		gcount=gcount+1
		!calculate the distance measure
    		do i=1,k
    			if(.not.doCal(i)) cycle
			do j=1,npt
				mdis(i,j)=MahaDis(ndim,pt(:,j),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
				h(i,j)=volx(i)*mdis(i,j)/nptx(i)
			enddo
    		enddo
    		
!		do x=1,npt
!			i=clusterx(x) !cluster in which point x lies
!			
!			if(doCal(i)) mdis(x,i)=MahaDis(ndim,pt(:,x),meanx(i,:),invcovx(i,:,:),kfacx(i)*effx(i))
!			d2=kfacx(i)*effx(i)
!			d1=(((nptx(i)/(nptx(i)-1.d0))**3)*(1.d0-mdis(x,i)/(nptx(i)-1.d0))-1.d0)*sqrt(detcovx(i)*(d2**(ndim)))/ &
!    				(((nptx(i)/(nptx(i)-1.d0))**1.5d0)*sqrt(1.d0-mdis(x,i)/(nptx(i)-1.d0))+1.d0)
!			
!			do j=1,k
!				if(j==i) then
!					h(j,x)=0.d0
!				else
!					if(doCal(j)) mdis(x,j)=MahaDis(ndim,pt(:,x),meanx(j,:),invcovx(j,:,:),kfacx(j)*effx(j))
!					d2=kfacx(j)*effx(j)
!					h(j,x)=(((nptx(j)/(nptx(j)+1.d0))**3)*(1.d0+mdis(x,j)/(nptx(j)+1.d0))-1.d0)*sqrt(detcovx(j)*(d2**(ndim)))/ &
!    						(((nptx(j)/(nptx(j)+1.d0))**1.5d0)*(sqrt(1.d0+mdis(x,j)/(nptx(j)+1.d0)))+1.d0)
!					h(j,x)=(h(j,x)+d1)!/dble(nptx(i)+nptx(j))
!				endif
!			enddo
!    		enddo
    
    		!now assign point j to the ellipsoid i such that h(i,j) is min
    		flag=.false.
    		do j=1,npt
    			indx=minloc(h(1:k,j))
			if(clusterx(j)/=indx(1)) then
				doCal(clusterx(j))=.true.
				doCal(indx(1))=.true.
				flag=.true.
				nptx(clusterx(j))=nptx(clusterx(j))-1
				nptx(indx(1))=nptx(indx(1))+1
				clusterx(j)=indx(1)
			endif
    		enddo
		
		!if a cluster has less than min_pt points then kill it
		if(flag) then
			do j=1,k
				if(nptx(j)<min_pt) then
					do i1=1,npt
						if(clusterx(i1)==j) then
							h(j,i1)=1.d99
							indx=minloc(h(1:k,i1))
							doCal(indx(1))=.true.
							nptx(j)=nptx(j)-1
							nptx(indx(1))=nptx(indx(1))+1
							clusterx(i1)=indx(1)
						endif
					enddo
					doCal(j)=.false.
					mdis(j,1:npt)=1.d99
					volx(j)=0.d0
					kfacx(j)=0.d0
					evalx(j,:)=0.d0
				endif
			enddo
			exit
		endif

		if(.not.flag) then
			!no point reassigned
			!check if clustering resulted in the new cluster having more than min_pt points
			if(nptx(k)>0 .and. nptx(k)<npt) then
				!check the boundary point of each ellipsoid
				do i=1,k
					if(nptx(i)==min_pt) cycle
			
					!find the boundary point
					d1=0.d0
					do j=1,npt
						if(clusterx(j)==i) then
							if(mdis(i,j)>d1) then
								d1=mdis(i,j)
								i1=j
							endif
						endif
					enddo
			
					!now calculate delF for all the other ellipsoids
					do j=1,k
						if(i==j .or. nptx(j)==0) cycle
						d1=kfacx(i)*effx(i)
						d2=kfacx(j)*effx(j)
						if(delF(ndim,nptx(i),nptx(j),mdis(i,i1),mdis(j,i1),detcovx(i), &
						detcovx(j),d1,d2)<0.) then
							doCal(i)=.true.
							doCal(j)=.true.
							flag=.true.
							nptx(i)=nptx(i)-1
							nptx(j)=nptx(j)+1
							clusterx(i1)=j
							exit
						endif
					enddo
					if(flag) exit
				enddo
			endif
		endif
		
		if(flag) then
			!point reassigned, recalculate the ellipsoid properties
			count=count+1
			!if no better partition found in the past 5 iterations then return
			if(count==20) exit
			if(nptx(k)>0 .and. nptx(k)<npt) then
				do i=1,k
					if(.not.doCal(i)) cycle
					i1=0
					do j=1,npt
						if(clusterx(j)==i) then
							i1=i1+1
							ptk(1:ndim,i1)=pt(1:ndim,j)
						endif
					enddo
					d1=pVol*dble(i1)/dble(npt)
					call CalcEllProp(i1,ndim,ptk(1:ndim,1:i1),meanx(i,:),covmatx(i,:,:), &
						invcovx(i,:,:),tmatx(i,:,:),evecx(i,:,:),evalx(i,:),detcovx(i),kfacx(i), &
						effx(i),volx(i),d1,.false.)
				enddo
				!check if this resulted in a better partition
				if(sum(volx(1:k))<sum(volk(1:k))) then
					!reset the counter
					count=0
					!save the best partition
					nptk=nptx
    					cluster=clusterx
    					meank=meanx
    					covmatk=covmatx
	    				invcovk=invcovx
    					tmatk=tmatx
    					evalk=evalx
    					eveck=evecx
    					kfack=kfacx
	    				effk=effx
    					detcovk=detcovx
    					volk=volx
				
					if(sum(volx(1:k))<fVol) then
						Dmeans=0
						if(abs(sum(volx(1:k))-pVol)/pVol<0.01) return
					else
						Dmeans=1
					endif
				endif
			endif
		else
			exit
		endif
	enddo
	if(Dmeans==0) exit
    enddo
    
    deallocate( h, mdis, ptk, mu_tmp, meanx, covmatx, invcovx, tmatx, evalx, evecx, kfacx, effx, detcovx, volx )
    deallocate( doCal, broken )
    deallocate( nptx, clusterx, cluster2 )
    
	
  end function Dmeans  
 
!----------------------------------------------------------------------
  !calculate the variation in weighted average of e-tightness functions of 2 
  !components for re-assigning a point from component 1 to 2
  !Choi, Wang & Kim, Eurographics 2007, vol. 26, No. 3
  function delF(ndim,n1,n2,mdis1,mdis2,detcov1,detcov2,kfac1,kfac2)
    implicit none
    
    !input variables
    integer ndim !dinemsionality
    integer n1,n2 !no. of points in each ellipsoid
    double precision mdis1,mdis2 !Mahalanobis distance of the point from both ellipsoids
    double precision detcov1,detcov2 !determinant of the covariance matrices of the ellipsoids
    double precision kfac1,kfac2 !overall enlargement factors of the ellipsoids
    
    !output variables
    double precision delF
    
    
    delF=(((n1/(n1-1.))**3)*(1.-mdis1/(n1-1.))-1.)*sqrt(detcov1*(kfac1**(ndim)))/ &
    	(((n1/(n1-1.))**1.5)*sqrt(1.-mdis1/(n1-1.))+1.)
    delF=delF+(((n2/(n2+1.))**3)*(1.+mdis2/(n2+1.))-1.)*sqrt(detcov2*(kfac2**(ndim)))/ &
    	(((n2/(n2+1.))**1.5)*(sqrt(1.+mdis2/(n2+1.)))+1.)
	
    delF=delF/(dble(n1+n2))
	
  end function delF  
 
!----------------------------------------------------------------------

end module kmeans_clstr
