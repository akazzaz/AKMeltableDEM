#ifndef PROFILE_H
#define PROFILE_H
//

//

class Cslice
{
	public:
	//AK mod start - Volume fractions
	double melt_frac;	/** Melt volume fraction in slice.*/
	double solid_frac;	/** Solid volume fraction in slice.*/
	double void_frac;	/** Void volume fraction in slice.*/
	//AK mod end - Volume fractions
	//AK mod start - Stress matrix
	Cmatrix stress;
	//AK mod end - Stress matrix
	Cvector V;						/**<Average particle velocity.*/
	double frac_solid; 				/**<Average solid fraction.*/
	double T;						/**<Average particle temperature.*/
	double Y; 						/**< y position of the middle of the slice.*/
	double vol_slice;				/**< Volume of the slice.*/
	std::vector <Cparticle *> P;			/**< list of particle (pointers) that have a piece of volume in the slice.*/
	std::vector <double> volume_in_slice; /**< list of particle volume that belong to the slice.*/
	void average();
	Cslice(){frac_solid=0;T=0;};//initialise in constructor	
	friend ofstream &operator<<(ofstream &,Cparticle);		/**< Print the particle's data into a file.*/ 
};



class Cprofile  
{
	
public:
	Cconfig config;
	Cprofile(double step_target, Cconfig conf);
	double step;
	std::vector <Cslice> slice;
};
#endif
