#ifndef INC_REQ_newGA
#define INC_REQ_newGA
#include "newGA.hh"
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <iostream>
	#include <regex>
	#include <boost/regex.hpp>

using namespace std;

skeleton newGA
{

	// Problem ---------------------------------------------------------------
	Problem::Problem ():_dimension(0),_cantidad_vehiculos(NULL),_largo_mapa(NULL),_ancho_mapa(NULL),_autonomia_vehiculo(NULL),_mapa(NULL)
	{}


	ostream& operator<< (ostream& os, const Problem& pbm)
	{
		os << endl << endl << "Number of Variables " << pbm._dimension
		   << endl;
		return os;
	}

	istream& operator>> (istream& is, Problem& pbm)
	{
		return is;
	}

	int Problem::loadInstance(Problem& pbm, char* ruta) const{
		char buffer[MAX_BUFFER];
		string line, subline;
		ifstream readFile(ruta);
		//cantidad vehiculos
		getline(readFile,line	);
		stringstream iss(line);
		getline(iss,subline, ' ');
		int n = subline.length();
		char char_array[n+1];
		strcpy(char_array, subline.c_str());
		pbm._dimension = atoi(char_array);
		pbm._cantidad_vehiculos = pbm._dimension;
		pbm._autonomia_vehiculo = new int[pbm._cantidad_vehiculos];
		pbm._movimientos_restantes_vehiculo = new int[pbm._cantidad_vehiculos];
		pbm._total_explorado = 0;
		// mapa _largo
		getline(readFile,line	);
		stringstream iss2(line);
		getline(iss2,subline, ' ');
		n = subline.length();
		char_array[n+1];
		strcpy(char_array, subline.c_str());
		pbm._largo_mapa = atoi(char_array);
		//mapa ancho
		getline(readFile,line	);
		stringstream iss3(line);
		getline(iss3,subline, ' ');
		n = subline.length();
		char_array[n+1];
		strcpy(char_array, subline.c_str());
		pbm._ancho_mapa = atoi(char_array);
		pbm._mapa = new int*[pbm._largo_mapa];
		pbm._mapa_explorado = new int*[pbm._largo_mapa];
		for (int l = 0; l < pbm._largo_mapa; l++) {
			pbm._mapa[l] = new int[pbm._ancho_mapa];
			pbm._mapa_explorado[l] = new int[pbm._ancho_mapa];
		}
		//autonomia
		getline(readFile,line	);
		stringstream iss4(line);
		int iterVehiculo = 0;
		while (getline(iss4,subline, ' ')){
			n = subline.length();
			char_array[n+1];
			strcpy(char_array, subline.c_str());
			pbm._autonomia_vehiculo[iterVehiculo] = atoi(char_array);
			pbm._movimientos_restantes_vehiculo[iterVehiculo] = atoi(char_array);
			iterVehiculo ++;
		}
		//matriz
		int fila = 0;
		pbm._cantidad_zonas =0;
		while(getline(readFile,line	))   {
			stringstream iss5(line);
			int columna = 0;
			while (getline(iss5,subline, ' ')){
				n = subline.length();
				char_array[n+1];
				strcpy(char_array, subline.c_str());
				int tipo_zona = atoi(char_array);
				if (tipo_zona > 0){
					pbm._cantidad_zonas ++;
				}
				pbm._mapa[fila][columna] = tipo_zona;
				pbm._mapa_explorado[fila][columna] = 0;
				columna ++;
			}
			fila ++;
		}

		//log
		printf("%s", "dimension: ");
		printf("%d\n", pbm._dimension);
		printf("%s", "_cantidad_vehiculos: ");
		printf("%d\n", pbm._cantidad_vehiculos);
		printf("%s", "_largo_mapa: ");
		printf("%d\n", pbm._largo_mapa);
		printf("%s", "_ancho_mapa: ");
		printf("%d\n", pbm._ancho_mapa);
		printf("%s", "autonomia: ");
		for(int i = 0; i < pbm._cantidad_vehiculos; i++){
			printf("%d-", pbm._autonomia_vehiculo[i]);
		}
		printf("%s\n", "");
		printf("%s\n", "matriz: ");
		for(int j = 0; j < pbm._largo_mapa; j++){
			for(int k = 0; k < pbm._ancho_mapa; k++){
				printf("%d-", (pbm._mapa[j][k]));
			}
			printf("%s\n", "");
		}
		printf("%s\n", "");

		return 0;
	}


	bool Problem::operator== (const Problem& pbm) const
	{
		if (_dimension!=pbm.dimension()) return false;
		return true;
	}

	bool Problem::operator!= (const Problem& pbm) const
	{
		return !(*this == pbm);
	}

	Direction Problem::direction() const
	{
		 //return maximize;
		 return minimize;
	}

	int Problem::dimension() const
	{
		return _dimension;
	}

	int Problem::getCantidadVehiculos() const
	{
		return _cantidad_vehiculos;
	}


	int Problem::getLargoMapa() const
	{
		return _largo_mapa;
	}

	int Problem::getAnchoMapa() const
	{
		return _ancho_mapa;
	}

	int& Problem::getAutonomiaVehiculo(int index) const
	{
		return _autonomia_vehiculo[index];
	}

	int& Problem::getValorPosicionMapa(int largo,int ancho) const
	{
		return _mapa[largo][ancho];
	}

	int Problem::realizarMovimiento(int vehiculo)
	{
		_movimientos_restantes_vehiculo[vehiculo] -=1;
		return _movimientos_restantes_vehiculo[vehiculo];
	}

	int Problem::explorarZona(int largo,int ancho)
	{
		if (_mapa_explorado[largo][ancho] == 0){
			_total_explorado ++;
		}
		_mapa_explorado[largo][ancho] +=1;
		return _mapa_explorado[largo][ancho];
	}

	bool Problem::isLoadZone(int largo,int ancho) const
	{
		return (_mapa[largo][ancho] == 2);
	}

	bool Problem::isObstacle(int largo,int ancho) const
	{
		return (_mapa[largo][ancho] == 0);
	}

	bool Problem::isAlreadyExplored(int largo,int ancho) const
	{
		return (_mapa_explorado[largo][ancho] > 0);
	}

	bool Problem::objetivoCumplido() const
	{
		return (_cantidad_zonas / _total_explorado >= 90);
	}

	Problem::~Problem()
	{
	}
	// Solution --------------------------------------------------------------

	Solution::Solution (const Problem& pbm):_pbm(pbm),_var(pbm.dimension())
	{}

	const Problem& Solution::pbm() const
	{
		return _pbm;
	}

	Solution::Solution(const Solution& sol):_pbm(sol.pbm())
	{
		*this=sol;
	}

	istream& operator>> (istream& is, Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++)
			is >> sol._var[i];
		return is;
	}

	ostream& operator<< (ostream& os, const Solution& sol)
	{
		for (int i=0;i<sol.pbm().dimension();i++)
			os << " " << sol._var[i];
		return os;
	}

	NetStream& operator << (NetStream& ns, const Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns << sol._var[i];
		return ns;
	}

	NetStream& operator >> (NetStream& ns, Solution& sol)
	{
		for (int i=0;i<sol._var.size();i++)
			ns >> sol._var[i];
		return ns;
	}

 	Solution& Solution::operator= (const Solution &sol)
	{
		_var=sol._var;
		return *this;
	}

	bool Solution::operator== (const Solution& sol) const
	{
		if (sol.pbm() != _pbm) return false;
		for(int i = 0; i < _var.size(); i++)
			if(_var[i] != sol._var[i]) return false;
		return true;
	}

	bool Solution::operator!= (const Solution& sol) const
	{
		return !(*this == sol);
	}

	void Solution::initialize()
	{
		// double tiempo_empleado[_pbm.cantidadempleados()+1];
		// for (int k=1;k<=_pbm.cantidadempleados();k++){
		// 	tiempo_empleado[k] = 0.0;
		// }
		// for (int i=0;i<_pbm.dimension();i++){
		// 	bool encontre = false;
		// 	while (! encontre){
		// 		int indice = rand_int(1,_pbm.cantidadempleados());
		// 		double aux = tiempo_empleado[indice];
		// 		aux += _pbm.gettiempotarea(i)/(0.5 +_pbm.gethabilidad(indice-1));
		// 		if (isValid(aux, _pbm.getdisponibilidad(indice-1), _pbm.limiteproyecto())){
		// 				_var[i]= indice;
		// 				tiempo_empleado[indice] = aux;
		// 				encontre = true;
		// 		}
		// 	}
		//
		// }


	}

	// bool Solution::isValid(double coef, int disponibilidad_empleado, int limite_proyecto){
	// 	double cantidad_dias_empleado = coef/disponibilidad_empleado;
	// 	double limite_proyecto_double = limite_proyecto * 1.0;
	// 	return (limite_proyecto_double >= cantidad_dias_empleado);
	// 	// return (limite_proyecto >= (int) ceil(cantidad_dias_empleado));
	// }

	double Solution::fitness ()
	{
    // double fitness = 0.0;
		// double tiempo_total = 0.0;
		// int limite_proyecto = _pbm.limiteproyecto();
		// double tiempo_empleado[_pbm.cantidadempleados()+1];
		// for (int k=1;k<=_pbm.cantidadempleados();k++){
		// 	tiempo_empleado[k] = 0.0;
		// }
		// for (int i=0;i<_var.size();i++){//recorro tupla solucion(tareas)
		// 	int indice = _var[i]; //Indice de empleado
		// 	tiempo_empleado[indice] += _pbm.gettiempotarea(i)/(0.5 +_pbm.gethabilidad(indice-1));
		// }
		// for(int j=1;j<=_pbm.cantidadempleados();j++){
		// 	int disponibilidad_empleado = _pbm.getdisponibilidad(j-1);
		// 	int sueldo_empleado = _pbm.getsueldo(j-1);
		// 	double cantidad_dias_empleado = tiempo_empleado[j]/disponibilidad_empleado;
		// 	fitness +=(int) ceil(cantidad_dias_empleado) * sueldo_empleado;
		// }
		// return fitness;
	}

	char *Solution::to_String() const
	{
		return (char *)_var.get_first();
	}

	void Solution::to_Solution(char *_string_)
	{
		int *ptr=(int *)_string_;
		for (int i=0;i<_pbm.dimension();i++)
		{
			_var[i]=*ptr;
			ptr++;
		}
	}

	unsigned int Solution::size() const
	{
		return (_pbm.dimension() * sizeof(int));
	}


	int& Solution::var(const int index)
	{
		return _var[index];
	}


	Rarray<int>& Solution::array_var()
	{
		return _var;
	}

	Solution::~Solution()
	{}

	// UserStatistics -------------------------------------------------------

	UserStatistics::UserStatistics ()
	{}

	ostream& operator<< (ostream& os, const UserStatistics& userstat)
	{
		os << "\n---------------------------------------------------------------" << endl;
		os << "                   STATISTICS OF TRIALS                   	 " << endl;
		os << "------------------------------------------------------------------" << endl;

		for (int i=0;i< userstat.result_trials.size();i++)
		{
			os << endl
			   << userstat.result_trials[i].trial
			   << "\t" << userstat.result_trials[i].best_cost_trial
			   << "\t\t" << userstat.result_trials[i].worst_cost_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_evaluation_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].nb_iteration_best_found_trial
			   << "\t\t\t" << userstat.result_trials[i].time_best_found_trial
			   << "\t\t" << userstat.result_trials[i].time_spent_trial;
		}
		os << endl << "------------------------------------------------------------------" << endl;
		return os;
	}

	UserStatistics& UserStatistics::operator= (const UserStatistics& userstats)
	{
		result_trials=userstats.result_trials;
		return (*this);
	}

	void UserStatistics::update(const Solver& solver)
	{
		if( (solver.pid()!=0) || (solver.end_trial()!=true)
		  || ((solver.current_iteration()!=solver.setup().nb_evolution_steps())
		       && !terminateQ(solver.pbm(),solver,solver.setup())))
			return;

		struct user_stat *new_stat;

		if ((new_stat=(struct user_stat *)malloc(sizeof(struct user_stat)))==NULL)
			show_message(7);
		new_stat->trial         		 		 = solver.current_trial();
		new_stat->nb_evaluation_best_found_trial = solver.evaluations_best_found_in_trial();
		new_stat->nb_iteration_best_found_trial  = solver.iteration_best_found_in_trial();
		new_stat->worst_cost_trial     		 	 = solver.worst_cost_trial();
		new_stat->best_cost_trial     		 	 = solver.best_cost_trial();
		new_stat->time_best_found_trial		 	 = solver.time_best_found_trial();
		new_stat->time_spent_trial 		 		 = solver.time_spent_trial();

		result_trials.append(*new_stat);
	}

	void UserStatistics::clear()
	{
		result_trials.remove();
	}

	UserStatistics::~UserStatistics()
	{
		result_trials.remove();
	}

// Intra_operator  --------------------------------------------------------------

	Intra_Operator::Intra_Operator(const unsigned int _number_op):_number_operator(_number_op),probability(NULL)
	{}

	unsigned int Intra_Operator::number_operator() const
	{
		return _number_operator;
	}

	Intra_Operator *Intra_Operator::create(const unsigned int _number_op)
	{
		switch (_number_op)
		{
			case 0: return new Crossover;break;
			case 1: return new Mutation();break;
		}
	}

	ostream& operator<< (ostream& os, const Intra_Operator& intra)
	{
		switch (intra.number_operator())
		{
			case 0: os << (Crossover&)intra;break;
			case 1: os << (Mutation&)intra;break;
		}
		return os;
	}

	Intra_Operator::~Intra_Operator()
	{}

//  Crossover:Intra_operator -------------------------------------------------------------

	Crossover::Crossover():Intra_Operator(0)
	{
		probability = new float[1];
	}

	void Crossover::cross(Solution& sol1,Solution& sol2) const
	{
		// bool iguales = true;
		// for(int y=0;y <sol1.pbm().dimension(); y++){
		// 	if (sol1.var(y) != sol2.var(y)){
		// 		iguales = false;
		// 		break;
		// 	}
		// }
		// if(iguales){
		// 	sol2.initialize();
		// }else{
		// 	int i=0;
		// 	Rarray<int> aux(sol1.pbm().dimension());
		// 	Rarray<int> aux2(sol2.pbm().dimension());
		// 	aux=sol2.array_var();
		// 	aux2=sol1.array_var();
		// 	bool end = false;
		// 	double tiempo_empleado[sol1.pbm().cantidadempleados() +1];
		// 	double tiempo_empleado2[sol1.pbm().cantidadempleados() +1];
		// 	double coef =0;
		// 	for (int k=1;k<=sol1.pbm().cantidadempleados();k++){
		// 		tiempo_empleado[k] = 0.0;
		// 		tiempo_empleado2[k] = 0.0;
		// 	}
		// 	double tiempo =0;
		// 	int limit=rand_int((sol1.pbm().dimension()/2)+1,sol1.pbm().dimension()-1);
		// 	int limit2=rand_int(0,limit-1);
		// 	for (i=0;i<limit2;i++){
		// 		sol2.var(i)=sol1.var(i);
		// 		coef = sol2.pbm().gettiempotarea(i)/(0.5 +sol2.pbm().gethabilidad(sol2.var(i)-1));
		// 		tiempo_empleado2[sol2.var(i)] += coef;
		// 	}
		//
		// 	for (int i=0;i<limit2;i++){
		// 		sol1.var(i)=aux[i];
		// 		coef = sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(sol1.var(i)-1));
		// 		tiempo_empleado[sol1.var(i)] += coef;
		// 	}
		//
		// 	for (i= limit2; i<limit;i++){
		// 		tiempo = tiempo_empleado2[sol2.var(i)];
		// 		tiempo += sol2.pbm().gettiempotarea(i)/(0.5 +sol2.pbm().gethabilidad(sol2.var(i)-1));
		// 		if(sol2.isValid(tiempo, sol2.pbm().getdisponibilidad(sol2.var(i)-1), sol2.pbm().limiteproyecto())){
		// 			tiempo_empleado2[sol2.var(i)] = tiempo;
		// 		}else{
		// 			tiempo = tiempo_empleado2[sol1.var(i)];
		// 			tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(sol1.var(i)-1));
		// 			if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(sol1.var(i)-1), sol1.pbm().limiteproyecto())){
		// 				tiempo_empleado2[sol1.var(i)] = tiempo;
		// 				sol2.var(i) = sol1.var(i);
		// 			}
		// 			else{
		// 				end = true;
		// 				break;
		// 			}
		// 		}
		// 		if(!end){
		// 			//SOLUCION1
		// 			tiempo = tiempo_empleado[sol1.var(i)];
		// 			tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(sol1.var(i)-1));
		// 			if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(sol1.var(i)-1), sol1.pbm().limiteproyecto())){
		// 				tiempo_empleado[sol1.var(i)] = tiempo;
		// 			}else{
		// 				tiempo = tiempo_empleado[aux[i]];
		// 				tiempo += sol2.pbm().gettiempotarea(i)/(0.5 +sol2.pbm().gethabilidad(aux[i]-1));
		// 				if(sol2.isValid(tiempo, sol2.pbm().getdisponibilidad(aux[i]-1), sol2.pbm().limiteproyecto())){
		// 					tiempo_empleado[aux[i]] = tiempo;
		// 					sol1.var(i) = aux[i];
		// 				}else{
		// 					end = true;
		// 					break;
		// 				}
		// 			}
		// 		}
		// 	}
		// 	if(!end){
		// 		for (i=limit;i<sol1.pbm().dimension();i++){
		// 			tiempo = tiempo_empleado2[sol1.var(i)];
		// 			double _habilidad_empleado = sol1.pbm().gethabilidad(sol1.var(i)-1);
		// 			tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +_habilidad_empleado);
		// 			if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(sol1.var(i)-1), sol1.pbm().limiteproyecto())){
		// 				tiempo_empleado2[sol1.var(i)] = tiempo;
		// 				sol2.var(i) = sol1.var(i);
		// 			}else{
		// 				tiempo = tiempo_empleado2[sol2.var(i)];
		// 				tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(sol2.var(i)-1));
		// 				if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(sol2.var(i)-1), sol1.pbm().limiteproyecto())){
		// 					tiempo_empleado2[sol2.var(i)] = tiempo;
		// 				}else{
		// 					end = true;
		// 					break;
		// 				}
		// 			}
		// 		}
		// 	}
		// 	if(!end){
		// 		for (i=limit;i<sol1.pbm().dimension();i++){
		// 			tiempo = tiempo_empleado[aux[i]];
		// 			tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(aux[i]-1));
		// 			if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(aux[i]-1), sol1.pbm().limiteproyecto())){
		// 				sol1.var(i)=aux[i];
		// 				tiempo_empleado[sol1.var(i)] = tiempo;
		// 			}else{
		// 				tiempo = tiempo_empleado[sol1.var(i)];
		// 				tiempo += sol1.pbm().gettiempotarea(i)/(0.5 +sol1.pbm().gethabilidad(sol1.var(i)-1));
		// 				if(sol1.isValid(tiempo, sol1.pbm().getdisponibilidad(sol1.var(i)-1), sol1.pbm().limiteproyecto())){
		// 					tiempo_empleado[sol1.var(i)] = tiempo;
		// 				}else{
		// 					end = true;
		// 					break;
		// 				}
		//
		// 			}
		// 		}
		// 	}
		// 	if(end){ //no hago el cruzamiento
		// 		sol2.array_var()=aux;
		// 		sol1.array_var()=aux2;
		// 	}
		// }



	}

	void Crossover::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i+1<sols.size();i=i+2)
		 	if (rand01()<=probability[0]) cross(*sols[i],*sols[i+1]);
	}

	ostream& operator<< (ostream& os, const Crossover&  cross)
	{
		 os << "Crossover." << " Probability: "
                    << cross.probability[0]
		    << endl;
		 return os;
	}

	void Crossover::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_crossover_probability",(char *)probability,1,sizeof(float));
	}

	void Crossover::UpdateFromState(const StateCenter& _sc)
	{
		 unsigned long nbytes,length;
		 _sc.get_contents_state_variable("_crossover_probability",(char *)probability,nbytes,length);
	}

	void Crossover::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f ",&op,&probability[0]);
		assert(probability[0]>=0);
	}

	Crossover::~Crossover()
	{
		delete [] probability;
	}

	//  Mutation: Sub_operator -------------------------------------------------------------

	Mutation::Mutation():Intra_Operator(1)
	{
		probability = new float[2];
	}

	void Mutation::mutate(Solution& sol) const
	{
		if (rand01()<=probability[1]){
	        // int i = rand_int(0, sol.pbm().dimension() -1);
					// bool encontre = false;
					// int iter = 0;
					//
					// while (! encontre && iter < 20){
					// 	int indice = rand_int(1,sol.pbm().cantidadempleados());
					// 	if (indice != sol.var(i)){
					// 		double tiempo_empleado = 0;
					// 		for (int j =0; j < sol.pbm().dimension(); j++){
					// 			if(sol.var(j)== indice && j != i){
					// 				tiempo_empleado += sol.pbm().gettiempotarea(j)/(0.5 +sol.pbm().gethabilidad(indice-1));
					// 			}
					// 		}
					// 		tiempo_empleado += sol.pbm().gettiempotarea(i)/(0.5 +sol.pbm().gethabilidad(indice-1));
					// 		if (sol.isValid(tiempo_empleado, sol.pbm().getdisponibilidad(indice-1), sol.pbm().limiteproyecto())){
					// 				sol.var(i)= indice;
					// 				encontre = true;
					// 		}
					// 	}
					// 	iter ++;
					// }
	    }
	}

	void Mutation::execute(Rarray<Solution*>& sols) const
	{
		for (int i=0;i<sols.size();i++)
			if(rand01() <= probability[0])	mutate(*sols[i]);
	}

	ostream& operator<< (ostream& os, const Mutation&  mutation)
	{
		os << "Mutation." << " Probability: " << mutation.probability[0]
		   << " Probability1: " << mutation.probability[1]
		   << endl;
		return os;
	}

	void Mutation::setup(char line[MAX_BUFFER])
	{
		int op;
		sscanf(line," %d %f %f ",&op,&probability[0],&probability[1]);
		assert(probability[0]>=0);
		assert(probability[1]>=0);
	}

	void Mutation::RefreshState(const StateCenter& _sc) const
	{
		_sc.set_contents_state_variable("_mutation_probability",(char *)probability,2,sizeof(probability));
	}

	void Mutation::UpdateFromState(const StateCenter& _sc)
	{
		unsigned long nbytes,length;
		_sc.get_contents_state_variable("_mutation_probability",(char *)probability,nbytes,length);
	}

	Mutation::~Mutation()
	{
		delete [] probability;
	}

// StopCondition_1 -------------------------------------------------------------------------------------

	StopCondition_1::StopCondition_1():StopCondition()
	{}

	bool StopCondition_1::EvaluateCondition(const Problem& pbm,const Solver& solver,const SetUpParams& setup)
	{
		bool fin = ((int)solver.current_iteration() == 10000);
		return fin;
	}

	StopCondition_1::~StopCondition_1()
	{}

	//------------------------------------------------------------------------
	// Specific methods ------------------------------------------------------
	//------------------------------------------------------------------------

	bool terminateQ (const Problem& pbm, const Solver& solver,
			 const SetUpParams& setup)
	{
		StopCondition_1 stop;
		return stop.EvaluateCondition(pbm,solver,setup);
	}
}
#endif
