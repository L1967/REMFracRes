
#include <primitive_mesh/contrainte_geo.h>
#include <iomanip>
#include <fstream>
#include <memory>
#include <random>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <cstdlib>


namespace fractures_intersect {

//-------------------------------------------------------
// Version 1.0 modifiee - Pascal SIEGEL - 2000
//-------------------------------------------------------

const long MAXCAR = 120;
//-------------------------------------------------------
// Classe Ccontrainte_Geo :
// regroupe les contraintes geometriques et geologiques
// discretis�es ou non 
//
// compos� d'un tableau de polyedres
// + tableau de  plan 
// + tableau de primitives
//
// Version 1.0 Pascal SIEGEL - 2000
//-------------------------------------------------------

//-------------------------------------------
// constructeur/destructeur de la classe Ccontrainte_Geo
//-------------------------------------------
Ccontrainte_Geo::Ccontrainte_Geo() 
{
}

Ccontrainte_Geo::~Ccontrainte_Geo()
{
}

//-----------------------------------------------------
// Attributs
//-----------------------------------------------------
void Ccontrainte_Geo::Set_domaine(TBoite dom) {	domaine = dom; }
TBoite	Ccontrainte_Geo::Get_domaine() const { return domaine; }

//-------------------------------------------
// Elimination d'une primitive du tableau
// des primitives
//-------------------------------------------
void Ccontrainte_Geo::Del_primitives(int i)
{
	Del_primitive(tab_primitives,i);
}

//-------------------------------------------
// Elimination d'un polyedre du tableau
// des polyedres
//-------------------------------------------
void Ccontrainte_Geo::Del_polyedres(int i)
{
	Del_polyedre(tab_polyedres,i);
}

//-------------------------------------------
// Elimination d'un plan du tableau
// des plans
//-------------------------------------------
void Ccontrainte_Geo::Del_plans(int i)
{
	Del_plan(tab_plans,i);
}

//-------------------------------------------
// Operations relatives aux intersections
//-------------------------------------------
void Ccontrainte_Geo::intersect_all()
{
	Del_intersect();

	int	npoly = tab_polyedres.size(),
		nplan = tab_plans.size(),
		nprim = tab_primitives.size();

	if (!npoly && !nplan && !nprim) return;


	for (int i = 0 ; i < npoly ; i++)
	{
		
		for (int j = i+1 ; j < npoly ; j++)
		{
			tab_polyedres[i]->intersect(tab_polyedres[j].get());
		}

		for (int j = 0 ; j < nplan ; j++)
		{
			tab_polyedres[i]->intersect(tab_plans[j].get());
		}

		for (int j = 0 ; j < nprim ; j++)
		{
			tab_polyedres[i]->intersect(tab_primitives[j].get());
		}
		
	}

	for (int i = 0 ; i < nplan ; i++)
	{
		
		for (int j = i+1 ; j < nplan ; j++)
		{
			tab_plans[i]->intersect(tab_plans[j].get());
		}

		for (int j = 0 ; j < nprim ; j++)
		{
			tab_plans[i]->intersect(tab_primitives[j].get());
		}
		
	}

	Cprimitive	*primi,*primj;
	Tprimitive	ptypi,ptypj;

	for (int i = 0 ; i < nprim ; i++)
	{
		primi = tab_primitives[i].get();
		ptypi = primi->Get_type();
		for (int j = i+1 ; j < nprim ; j++)
		{
				
			primj = tab_primitives[j].get();
			ptypj = primj->Get_type();
			// intersection d'un polygone avec une polyligne
			// et non le contraire
			if ((ptypi == Polyligne && ptypj == Polygone)
				|| (ptypi == Point && ptypj != Point))
				primj->intersect(primi);
			else
				primi->intersect(primj);
			
		}
	}
	
	// debut debug info intersection
	Fusionne_primitives();
	scinde_lignes();
	// fin debug info intersection
}

void Ccontrainte_Geo::add_fracture_set(std::shared_ptr<fracture_set_geo> fracture_set) {


	tab_fractures_set_.push_back(fracture_set);

}


//-------------------------------------------
// Operations relatives aux intersections
//-------------------------------------------
void Ccontrainte_Geo::intersect_all_fracture_set(int number_of_thread)
{

	std::uniform_real_distribution<double> distribution(0.0,1.0);



	for (size_t si = 0 ; si < tab_fractures_set_.size() ; si++)
	{

		std::mt19937 generator(tab_fractures_set_[si]->seed_);

//		#pragma  omp  parallel for schedule(dynamic,5)  num_threads(number_of_thread) if (tab_fractures_set_[si]->tab_fractures_.size() > 10)

		////std::cout << "si = " << si << std::endl;
		for (size_t fi  = 0 ; fi < tab_fractures_set_[si]->tab_fractures_.size() ; fi++)
		{
			Cfaille* fracture_ref  = tab_fractures_set_[si]->tab_fractures_[fi].get();
			//std::cout << "fi = " << fi << std::endl;
			for (size_t sj = 0 ; sj < si ; sj++)
			{
				bool random = false;
				//std::cout << "sj = " << sj << " tab_fractures_set_[sj]->set_type_ " << tab_fractures_set_[sj]->set_type_ << " size " << tab_fractures_set_[si]->tab_type_intersection_.size() <<std::endl;
				if (tab_fractures_set_[si]->tab_type_intersection_[tab_fractures_set_[sj]->set_type_] == CROSSING) continue;
				if (tab_fractures_set_[si]->tab_type_intersection_[tab_fractures_set_[sj]->set_type_] == RANDOM) {
					random = true;
				}

				for (size_t fj = 0 ; fj < tab_fractures_set_[sj]->tab_fractures_.size() ; fj++)
				{
					Cfaille* fracture= tab_fractures_set_[sj]->tab_fractures_[fj].get();
					double prob = 1.0;
					if (random) prob = distribution(generator);

					if (prob > 0.5) {
						fracture_ref->intersect_plan_limitated(fracture,pt_tol);
					}
				}
			}

		}

	}
}



	
//---------------------------------------------
// On regarde dans les tableaux de primitives
// des differents objets si des polylignes
// sont discontinues ou non
//---------------------------------------------
void Ccontrainte_Geo::scinde_lignes()
{
	for (size_t i = 0 ; i < tab_polyedres.size() ; i++)
	{
		tab_polyedres[i]->scinde_lignes();
	}

	for (size_t i = 0 ; i < tab_plans.size() ; i++)
	{
		tab_plans[i]->scinde_lignes();
	}
	
	scinde_ligne(tab_primitives);
}

//-------------------------------------------
// Fusionne les primitives identiques
// issues d'intersection
//-------------------------------------------
void Ccontrainte_Geo::Fusionne_primitives()
{
	int	npoly = tab_polyedres.size();
	int nplan = tab_plans.size();


	for (int i = 0 ; i < npoly ; i++)
	{
		tab_polyedres[i]->Fusionne_primitives();
	}

	for (int i = 0 ; i < nplan ; i++)
	{
		tab_plans[i]->Fusionne_primitives();
	}
	
	Fusionne(tab_primitives);
}

//-------------------------------------------
// Elimine toutes les intersections trouv�es 
// auparavant
//-------------------------------------------
void Ccontrainte_Geo::Del_intersect()
{
	int	npoly = tab_polyedres.size(),
		nplan = tab_plans.size(),
		nprim = tab_primitives.size();

	if (!npoly && !nplan && !nprim) return;
	
	for (int i = npoly-1 ; i >= 0 ; i--)
	{
		if (!tab_polyedres[i]->tab_plan.size())
			Del_polyedres(i);

		else if (tab_polyedres[i] != NULL)
			tab_polyedres[i]->Del_intersect();
	}

	for (int i = nplan-1 ; i >= 0 ; i--)
	{
		if (!tab_plans[i]->tab_primitives.size())
			Del_plans(i);

		else if (tab_plans[i] != NULL)
			tab_plans[i]->Del_intersect();
	}

	for (int i = nprim-1 ; i >= 0 ; i--)
	{
		if (!tab_primitives.size())
			Del_primitive(tab_primitives,i);

		else if (tab_primitives[i] != NULL)
			tab_primitives[i]->Del_intersect();
	}
}


//-------------------------------------------
// Verifie que le nom donn� n'est pas celui
// d'une primitive existante
//-------------------------------------------

std::shared_ptr<Cbase_primitive> Ccontrainte_Geo::Recherche_par_nom(std::string nom)
{
	int	npoly = tab_polyedres.size(),
		nplan = tab_plans.size(),
		nprim = tab_primitives.size();
	
	for (int i = 0 ; i < npoly ; i++)
	{
		if (nom == tab_polyedres[i]->Get_nom()) 
			return tab_polyedres[i];

		for (size_t j = 0; j < tab_polyedres[i]->tab_plan.size() ; j++)
		{		
			if (nom == tab_polyedres[i]->tab_plan[j]->Get_nom()) 
				return tab_polyedres[i]->tab_plan[j];
		}
	}

	for (int i = 0 ; i < nplan ; i++)
	{
		if (nom == tab_plans[i]->Get_nom()) 
			return tab_plans[i];
	}

	for (int i = 0 ; i < nprim ; i++)
	{
		if (nom == tab_primitives[i]->Get_nom())
			return tab_primitives[i];
	}

	return NULL;
}


//-------------------------------------------
// d�truit la primitive correspondant au nom 
//-------------------------------------------

void Ccontrainte_Geo::Remove_by_name(std::string nom)
{

	int	npoly = tab_polyedres.size(),
		nplan = tab_plans.size(),
		nprim = tab_primitives.size();
	
	for (int i = 0 ; i < npoly ; i++)
	{
		if (nom==tab_polyedres[i]->Get_nom())
		{
			tab_polyedres.erase(tab_polyedres.begin()+i);
			return;
		}
	}		

	for (int i = 0 ; i < nplan ; i++)
	{
		if (nom==tab_plans[i]->Get_nom())
		{
			tab_plans.erase(tab_plans.begin()+i);
			return;
		}
	}

	for (int i = 0 ; i < nprim ; i++)
	{
		if (nom==tab_primitives[i]->Get_nom())
		{	
			Del_primitives(i);
			return ;
		}
	}

}


//-------------------------------------------
// d�truit la primitive correspondant au nom
//-------------------------------------------

void Ccontrainte_Geo::Remove_all()
{

	int	npoly = tab_polyedres.size(),
		nplan = tab_plans.size(),
		nprim = tab_primitives.size();

	for (int i = 0 ; i < npoly ; i++)
	{
		tab_polyedres.erase(tab_polyedres.begin()+i);

	}

	for (int i = 0 ; i < nplan ; i++)
	{

			tab_plans.erase(tab_plans.begin()+i);

	}

	for (int i = 0 ; i < nprim ; i++)
	{

			Del_primitives(i);
			return ;

	}
	for (size_t si = 0 ; si < tab_fractures_set_.size() ; si++)
	{
		tab_fractures_set_.erase(tab_fractures_set_.begin()+si);
	}
}


//-------------------------------------------
// ajoute une primitive au geo_contraintes
//-------------------------------------------

void Ccontrainte_Geo::Add_primitive(std::shared_ptr<Cbase_primitive> prim)
{
	Tprimitive type;

	type=prim->Get_type();


	if (type==Polyedre || type==Cylindre || type==Boite_Rect)
	{
		tab_polyedres.push_back(std::dynamic_pointer_cast<Cpolyedre> (prim));
		return;
	}

	if (type==Plan || type==Faille)
	{
		tab_plans.push_back(std::dynamic_pointer_cast<Cplan>(prim));
		return;
	}
		
	if (type==Point || type==Polyligne || type==Polygone)
	{
		tab_primitives.push_back(std::dynamic_pointer_cast<Cprimitive>(prim) );
		return;
	}
	
}




//-------------------------------------------
// change le polyedre contour
// l'ancien contour est remis dans le tableau
//-------------------------------------------


void Ccontrainte_Geo::Change_contour(std::shared_ptr<Cpolyedre> poly)
{
	std::shared_ptr<Cpolyedre> poly1;
	
	if (tab_polyedres.size()>0)
	{
		poly1=tab_polyedres[0];
		tab_polyedres[0]=poly;
		tab_polyedres.push_back(poly1);
	}
	else
	{
		tab_polyedres.push_back(poly);
	}
}

//-------------------------------------------
// operation sur le polyedre contour
//-------------------------------------------

void Ccontrainte_Geo::Set_countour(std::shared_ptr<Cpolyedre> poly)
{	


	if (tab_polyedres.size()>0)
	{
		tab_polyedres[0]=poly;
	}
	else
	{
		tab_polyedres.push_back(poly);
	}
}
	

std::shared_ptr<Cpolyedre> Ccontrainte_Geo::Get_countour()
{
	if (tab_polyedres.size()>0)
	{
		return tab_polyedres[0];
	}
	
	return std::shared_ptr<Cpolyedre>(nullptr);
}




void Ccontrainte_Geo::get_fracture_set_surface_mesh(int fracture_set_index, std::vector<Point3D>& tab_point, std::vector<std::vector<int>>& tab_triangle,
		std::vector<property<float>>& tab_triangle_fracture_property)
{
	if (fracture_set_index < tab_fractures_set_.size())  tab_fractures_set_[fracture_set_index]->get_fracture_set_surface_mesh(tab_point, tab_triangle, tab_triangle_fracture_property);
}



//-----------------------------------------------------
// operateur d'ecriture pour la classe Ccontrainte_Geo
//-----------------------------------------------------
std::ostream& operator<<(std::ostream& os, Ccontrainte_Geo* geo_imp)
{   

	os<<std::string("======================================================")<<"\n";
    os<<std::string("Fichier des contraintes geometriques et geologiques  :")<<"\n";
    os<<std::string("======================================================")<<"\n";
 
	os<<"\n";

	os<<std::string("Domaine d'etude : ")<<"\n";
	os<<geo_imp->domaine;
	os<<"\n";

	os<<std::string("Nombre de polyedres : ")<<geo_imp->tab_polyedres.size()<<"\n";
	os<<std::string("Nombre de plans ou failles : ")<<geo_imp->tab_plans.size()<<"\n";
	os<<std::string("Nombre de primitives autres : ")<<geo_imp->tab_primitives.size()<<"\n";

	if (geo_imp->tab_polyedres.size()>0)
	{
		os<<"\n";

		for(size_t i=0;i<geo_imp->tab_polyedres.size();i++) // ecriture des polyedres
		{
			os<<geo_imp->tab_polyedres[i];
		}
	}

	if (geo_imp->tab_plans.size()>0)
	{
		os<<"\n";

		for(size_t i=0;i<geo_imp->tab_plans.size();i++) // ecriture des plans et failles
		{ 
			os<<geo_imp->tab_plans[i];
		}
	}

	
	if (geo_imp->tab_primitives.size()>0)
	{
		os<<"\n";

		for(size_t i=0;i<geo_imp->tab_primitives.size();i++) // ecriture des primitives
		{
			os<<geo_imp->tab_primitives[i];
		}
	}

	return os;
} 

//-----------------------------------------------------
// operateur de lecture pour la classe Ccontrainte_Geo
//-----------------------------------------------------
std::istream& operator>>(std::istream& is, Ccontrainte_Geo* geo_imp)
{
/*	char d[MAXCAR];
	Tprimitive ptype;
	int nobjet;

	is.getline(d,MAXCAR,'\n');
	is.getline(d,MAXCAR,'\n');
	is.getline(d,MAXCAR,'\n');
	is.getline(d,MAXCAR,'\n');

	is.getline(d,MAXCAR,'\n');
	is>>geo_imp->domaine;
	is.getline(d,MAXCAR,'\n');
	

	is.getline(d,MAXCAR,':');
	is>>nobjet;
	is.getline(d,MAXCAR,'\n');
	geo_imp->tab_polyedres.resize(nobjet);

	is.getline(d,MAXCAR,':');
	is>>nobjet;
	is.getline(d,MAXCAR,'\n');
	geo_imp->tab_plans.resize(nobjet);

	is.getline(d,MAXCAR,':');
	is>>nobjet;
	is.getline(d,MAXCAR,'\n');
	geo_imp->tab_primitives.resize(nobjet);

	if (geo_imp->tab_polyedres.size()>0)  // lectures des polyedres
	{
		is.getline(d,MAXCAR,'\n');

		for(int i=0;i<geo_imp->tab_polyedres.size();i++) 
		{
			is>>ptype;

			else if (ptype == Boite_Rect)
				geo_imp->tab_polyedres[i]=new Cboite_rect;

			else
				geo_imp->tab_polyedres[i]=new Cpolyedre;

			is>>geo_imp->tab_polyedres[i];
		}
	}

	if (geo_imp->tab_plans.size()>0)   // lecture des plans et failles
	{
		is.getline(d,MAXCAR,'\n');

		for(int i=0;i<geo_imp->tab_plans.size();i++) 
		{	
			is>>ptype;

			if (ptype == Faille)
				geo_imp->tab_plans[i] = new Cfaille;

			else
				geo_imp->tab_plans[i] = new Cplan;
			
			is >> geo_imp->tab_plans[i];	
		}
	}

	if (geo_imp->tab_primitives.size()>0)	// lecture des primitives 3D
	{
		is.getline(d,MAXCAR,'\n');

		for(int i=0;i<geo_imp->tab_primitives.size();i++) 
		{
			is>>ptype;
	
			switch (ptype)
			{
				case Polygone: 

					geo_imp->tab_primitives[i]=new Cpolygone;
				
				break;

				case Polyligne:

					geo_imp->tab_primitives[i]=new Cpolyligne;

				break;

				case Point: 
				
					geo_imp->tab_primitives[i] = new Cprimitive;
					
				break;
			}

			is>>geo_imp->tab_primitives[i];	
		}
	}
	*/
	return is;
}


}


