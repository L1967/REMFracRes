#include <primitive_mesh/plan.h>
#include <primitive_mesh/primitive.h>
#include <vector>
#include <fstream>
#include <iostream>

#if !defined(LaFaille)
	#define LaFaille

namespace fractures_intersect {

//-------------------------------------------
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------

//-------------------------------------------
// classe Cfaille : public Cplan
// 
// equation du plan + limites externes
// et internes : contour + tableau de 
// primitives du plan 
// + pour la faille : azimut et pendage
//-------------------------------------------
struct TBoite;
class Cfaille : public Cplan
{
	// variables
public:

	double	azimut, // azimut : angle compris entre 0 et 360�
			pendage,// pendage : angle compris entre 0 et 90�
			ext2L_selon_lpgp,// extension 2L // lpgp
			ext2l_orthog_lpgp;// extension 2l orthogonale a lpgp

	int fracture_set_index_;
	double aperture_;

	double stress_zone_max_;

	bool	limite;	// plan limite par l'utilisateur sinon limite=intersection avec 
					// polyedre de contour 
	int	connection_relay_;
	double ratio_perm_;
	double slope_azimuth_;
	int familly_index_;

public:

	// constructeur, destructeur
	Cfaille(int fracture_set_index = 0.0 , Point3D pt_ref = Point3D(),double azim = 0.,
			double pend = 0.,double ext_selon_lpgp = 0.,
			double ext_ortho_lpgp = 0., double aperture = 0.0, double stress_zone_max = 0.0);

	// constructeur, destructeur
	Cfaille(Cpolygone* polygone);
	Cfaille(int fracture_set_index, double aperture, std::vector<Point3D>& poly, std::string name);
	
	Cfaille(const Cfaille& f);
	~Cfaille();
	
	// fonctions
	double	Get_azimut() const;
	double	Get_pendage() const;
	bool	Get_limite() const;
	double	Get_ext2L_selon_lpgp() const;
	double	Get_ext2l_orthog_lpgp() const;
	void	Set_azimut(double azim),
			Set_pendage(double pend),
			Set_limite(bool limit),
			Set_ext2L_selon_lpgp(double l),
			Set_ext2l_orthog_lpgp(double l);
	void	Set_SlopeAzimut(double azim);
	double	Get_SlopeAzimut() const;
	void	Set_RatioPermeability(double azim);
	double	Get_RatioPermeability() const;
	void	Set_aperture(double aperture);
	double	Get_aperture() const;
			
	// remet a jour les vecteurs qui definissent la faille
	virtual bool Actualise_repere(bool cal = false);
	
	bool	Set_faille_contour(TBoite *box);
	void	Creation_faille_finie();

	void simplify_to_box_poly();

    bool in_stress_zone(Point3D& pt);

    bool in_Plan(Point3D& pt);

	Cfaille& operator=(const Cfaille& f);

	void save(std::fstream& fout);
};

}

#endif
