//
ofstream & operator<<(ofstream &file,Cslice s)
{
	file<<s.Y<<"\t";
	// AK mod start - Volume fractions
	file<<s.melt_frac<<"\t";
	file<<s.solid_frac<<"\t";
	file<<s.void_frac<<"\t";
	// AK mod end - Volume fractions
	file<<s.frac_solid<<"\t";
	file<<s.V;
	file<<s.T<<"\t";
	// AK mod start - stress
	file<<s.stress.x[1][1]<<"\t";
	file<<s.stress.x[1][0]<<"\t";
	// AK mod start - stress
	file<<endl;	//new line
 	return file;
}

void Cslice::average(void)
{
	double volume_total_in_slice=0;
	// AK mod start - Calculate volume fractions
	double melt_vol=0;
	double solid_vol=0;
	double olapVs=0, olapVsi=0, olapVm=0;
	// AK mod end - Calculate volume fractions
	for(int ip=0;ip<P.size();ip++)
		{
			volume_total_in_slice+=volume_in_slice[ip];
			T+= P[ip]->T* volume_in_slice[ip];
			V+= (P[ip]->V* volume_in_slice[ip]);
			// AK mod start - Calculate volume fractions
			if (BRANCH != "CREATE"){
				solid_vol+=(P[ip]->RS/P[ip]->R)*volume_in_slice[ip];
				melt_vol+=(1-(P[ip]->RS/P[ip]->R))*volume_in_slice[ip];
			}
			else {solid_vol+=volume_in_slice[ip];}
			// AK mod end - Calculate volume fractions
		}
	frac_solid = volume_total_in_slice/vol_slice;
	T/=volume_total_in_slice;
	V/=volume_total_in_slice;
	// AK mod start - Calculate stress tensor
	Cmatrix sigma;
	for(int ic=0;ic<C.size();ic++) {
		olapVsi=PI/(12*C[ic]->dx)*pow((C[ic]->pA->RS+C[ic]->pB->RS-C[ic]->dx),2)*(pow(C[ic]->dx,2)+2*C[ic]->dx*C[ic]->pA->RS+2*C[ic]->dx*C[ic]->pB->RS-3*pow(C[ic]->pA->RS,2)-3*pow(C[ic]->pB->RS,2)+6*C[ic]->pA->RS*C[ic]->pB->RS);
     		olapVs+=olapVsi;
     		olapVm+=PI/(12*C[ic]->dx)*pow((C[ic]->pA->R+C[ic]->pB->R-C[ic]->dx),2)*(pow(C[ic]->dx,2)+2*C[ic]->dx*C[ic]->pA->R+2*C[ic]->dx*C[ic]->pB->R-3*pow(C[ic]->pA->R,2)-3*pow(C[ic]->pB->R,2)+6*C[ic]->pA->R*C[ic]->pB->R)-olapVsi;
		sigma+=(C[ic]->F| C[ic]->dX)*Cfrac[ic];
	}
	stress=sigma/volume_total_in_slice;
	// AK mod end - Calculate stress tensor
	// AK mod start - Calculate volume fractions
	melt_frac=((melt_vol-olapVm)>0.0?(melt_vol-olapVm):0.0)/vol_slice;
	solid_frac=((solid_vol-olapVs)>0.0?(solid_vol-olapVs):0.0)/vol_slice;
	void_frac=(1-melt_frac-solid_frac)>0.0?(1-melt_frac-solid_frac):0.0;
	// AK mod end - Calculate volume fractions
}




Cprofile::Cprofile(double step_target,Cconfig config) 
{
	double Nslice;
	
	//initiate the slices
	Nslice = (int)(config.cell.L.x[1]/step_target);
	step = config.cell.L.x[1]/Nslice; //trick to get same size of slice, close to step target
	
	Cslice temp;
	if(PSEUDO_2D) temp.vol_slice=step*config.cell.L.x[0]; //in that case, we assume that the z size is 1
	else  temp.vol_slice=step*config.cell.L.x[0]*config.cell.L.x[2];
	for(int i=0;i<Nslice;i++) 
	{
		temp.Y= i*step + step/2 - config.cell.L.x[1]/2; 
		slice.push_back(temp); //bluild the list of slice
	}
	//end initiate the slices

	//dispach_particle_in_slice();
	for(int ip=0; ip <config.P.size();ip++ )
		if(config.P[ip].AM_I_BOUNDARY==0)//only flowing particles
		{
			double y = config.P[ip].X.x[1] + config.cell.L.x[1]/2.;
			double R = config.P[ip].R;
	
			int imin = (int)((y-R)/step);
			int imax = (int)((y+R)/step);
		
			double a = -R;
			double b = (imin+1)*step-y;
			for(int islice=imin;islice<=imax;islice++)
				{
					int islice_periodic=islice;
					if(islice<0)					islice_periodic = islice+slice.size();
					else if(islice>=slice.size()) 	islice_periodic = islice-slice.size();
					
					if(islice_periodic<0|| islice_periodic>=slice.size())STOP("profile.cpp", "profile()","Bad slice index");
					
					slice[islice_periodic].P.push_back(&config.P[ip]);//get the particle into the slice
				
					if(islice==imax) b=R;	
					
					double vol = fabs(PI*(R*R*(b-a) - (b*b*b - a*a*a)/3.));
					
					slice[islice_periodic].volume_in_slice.push_back(vol);//get its volume belonging to the slice
					
					a=b;b+=step;
				}
				
		}
//end dispach_particle_in_slice();
// AK mod start - dispatch contacts in slice
	for(int ic=0; ic <config.C.size();ic++ )
		if(config.C[ic].pA->AM_I_BOUNDARY==0 && config.C[ic].pB->AM_I_BOUNDARY==0)//only flowing particles
		{
			double yA = config.C[ic].pA->X.x[1] + config.cell.L.x[1]/2.;
			double yB = config.C[ic].pB->X.x[1] + config.cell.L.x[1]/2.;
			
			double ymin = yA<yB?yA:yB;
			double ymax = yA<yB?yB:yA;
	
			int imin = (int)(ymin/step);
			int imax = (int)(ymax/step);
			
			for(int islice=imin;islice<=imax;islice++)
				{
					slice[islice].C.push_back(&config.C[ic]);//get the contact into the slice
				
					if (ymin>=islice*step && ymax>=(islice+1)*step){slice[islice].Cfrac.push_back(((islice+1)*step-ymin)/(ymax-ymin));}
					else if (ymin<=islice*step && ymax<=(islice+1)*step){slice[islice].Cfrac.push_back((ymin-islice*step)/(ymax-ymin));}
					else if (ymin<=islice*step && ymax>=(islice+1)*step){slice[islice].Cfrac.push_back(step/(ymax-ymin));}
					else {slice[islice].Cfrac.push_back(1.0);}
				}
				
		}
	//AK mod end - dispatch contact in slice

	for(int is=0;is<slice.size();is++) slice[is].average();
}
//
