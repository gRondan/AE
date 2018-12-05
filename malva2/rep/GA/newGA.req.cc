#ifndef INC_REQ_newGA
#define INC_REQ_newGA
#include "newGA.hh"
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <regex>
#include <boost/regex.hpp>
#include <tuple>
#include <vector>

using namespace std;

typedef vector< tuple<int,int> > tupleList;

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
		pbm._dimension = 10000;
		pbm._cantidad_vehiculos = atoi(char_array);
		pbm._autonomia_vehiculo = new int[atoi(char_array)];

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
		for (int l = 0; l < pbm._largo_mapa; l++) {
			pbm._mapa[l] = new int[pbm._ancho_mapa];
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
				columna ++;
			}
			fila ++;
		}

		//log
		// printf("%s", "dimension: ");
		// printf("%d\n", pbm._dimension);
		// printf("%s", "_cantidad_vehiculos: ");
		// printf("%d\n", pbm._cantidad_vehiculos);
		// printf("%s", "_largo_mapa: ");
		// printf("%d\n", pbm._largo_mapa);
		// printf("%s", "_ancho_mapa: ");
		// printf("%d\n", pbm._ancho_mapa);
		// printf("%s", "autonomia: ");
		// for(int i = 0; i < pbm._cantidad_vehiculos; i++){
		// 	printf("%d-", pbm._autonomia_vehiculo[i]);
		// }
		// printf("%s\n", "");
		// printf("%s\n", "matriz: ");
		// for(int j = 0; j < pbm._largo_mapa; j++){
		// 	for(int k = 0; k < pbm._ancho_mapa; k++){
		// 		printf("%d-", (pbm._mapa[j][k]));
		// 	}
		// 	printf("%s\n", "");
		// }
		// printf("%s\n", "");

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

	bool Problem::isLoadZone(int largo,int ancho) const
	{
		return (_mapa[largo][ancho] == 2 || isBasePosition(largo, ancho));
	}

	bool Problem::isObstacle(int largo,int ancho) const
	{
		return (_mapa[largo][ancho] == 0);
	}

	bool Problem::isExplorable(int largo, int ancho)const{
		return(largo >= 0 && ancho >= 0 && largo < _largo_mapa && ancho < _ancho_mapa && !isObstacle(largo, ancho));
	}

	int Problem::getCantidadZonas()const{
		return _cantidad_zonas;
	}

	bool Problem::isBasePosition(int largo, int ancho) const {
		return (largo == 0 && ancho == 0);
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

	vector<tupleList> Solution::getCaminos(){
		return caminos;
	}

	void Solution:: setCaminos(vector<tupleList> caminos2){
		caminos.swap(caminos2);
	}

	bool Solution::operator!= (const Solution& sol) const
	{
		return !(*this == sol);
	}

	int Solution:: total_explorado(){
		return _total_explorado;
	}

	void Solution::initialize()
	{
		//inicializo atributos de la solucion

		// printf("%s\n", "initialize");
		// printf("%d\n", _total_explorado);
		// printf("%d\n", caminos.size());

		std::vector<tupleList>().swap(caminos);
		_movimientos_restantes_vehiculo = new int[_pbm.getCantidadVehiculos()];
		for(int k = 0; k<_pbm.getCantidadVehiculos();k++){
			_movimientos_restantes_vehiculo[k] = _pbm.getAutonomiaVehiculo(k);
		}
		initializeMapaExplorado();

		//construyo los camino
		//la solucion testa compuesta con una lista que representa cada vehiculo donde cada elemento de la lista es una lista con las zonas del camino
		bool greedy = false;
		while(!objetivoCumplido()){
			for (int i = 0; i < _pbm.getCantidadVehiculos(); i ++){
				int totalExploradoInicial = _total_explorado;
				tupleList camino = construirCamino(i, greedy);

				// for(int j = 0; j < camino.size(); j++){
				// 	printf("%d-", get<0>(camino[j]));
				// 	printf("%d\n", get<1>(camino[j]));
				// }
				// printf("%s\n", "finCamino");
				if((camino.size() > 0) && (totalExploradoInicial < _total_explorado)){

					caminos.emplace_back(camino);
					std::vector<tuple<int,int>>().swap(camino);
					if(objetivoCumplido()){
						break;
					}
					greedy = false;
				}else{
					greedy = true;
				}

			}
		}
		convertCaminosToVar();
		std::vector<tupleList>().swap(caminos);
	}

	void Solution::initializeMapaExplorado(){
		_total_explorado = 0;
		_mapa_explorado = new int*[_pbm.getLargoMapa()];
		for (int l = 0; l < _pbm.getLargoMapa(); l++) {
			_mapa_explorado[l] = new int[_pbm.getAnchoMapa()];
		}
		for(int i = 0; i < _pbm.getLargoMapa(); i++){
			for(int j = 0; j < _pbm.getAnchoMapa(); j++){
				_mapa_explorado[i][j] = 0;
			}
		}
	}

	void Solution::convertCaminosToVar(){

			// printf("%s\n", "convertCaminosToVar");
		int iterVar = 0;
		bool caminoVacio = true;
			for (int i=0;i<caminos.size();i++){
				if(iterVar >= 9900){
					break;
				}
				tupleList camino = caminos[i];
				for(int j=0; j<camino.size(); j++){
					caminoVacio = false;
					if(iterVar >= 9900){
						break;
					}
					tuple<int,int> zona = camino[j];
					_var[iterVar] = get<0>(zona);
					iterVar ++;
					_var[iterVar] = get<1>(zona);6;
					iterVar ++;
				}
			//marco fin camino
			if(! caminoVacio){

				_var[iterVar] = -1;
				caminoVacio = true;
				iterVar ++;
			}
		}
		//fin de valores validos
		_var[iterVar] = -2;
		_var[iterVar +1] = -2;
		// for(int i = 0 ;i <_var.size(); i++){
		// 		printf("%d-", _var[i]);
		// }
		std::vector<tupleList>().swap(caminos);
		// printf("%s\n", "end convertCaminosToVar");
	}

	void Solution::convertVarToCaminos(){
		// printf("%s\n", "convertVarToCaminos");
		std::vector<tupleList>().swap(caminos);
		tupleList camino;
		bool caminoVacio = true;
		int inc = 1;
		for (int i=0;i<_var.size();i+=inc){
			if (_var[i] == -2){
				// fin de valores validos
				break;
			}
			if(_var[i] == -1){
				//fin del camino
				if(!caminoVacio){

					caminos.emplace_back(camino);
					std::vector<tuple<int,int>>().swap(camino);
				}
				caminoVacio = true;
				inc = 1;
			}else{
				caminoVacio = false;
				camino.emplace_back(_var[i], _var[i+1]);
				inc = 2;
			}

		}
		// printf("%s\n", "end convertVarToCaminos");
	}

	bool Solution::cumpleObjetivo(){
		// limpio el _mapa_explorado
		initializeMapaExplorado();
		//asigno los puntos explorados
		// printf("%s","caminos size " );
		// printf("%d\n", caminos.size());
		for(int i = 0; i<caminos.size(); i++){
			tupleList camino = caminos[i];

			for(int j = 0; j<camino.size(); j++){
				int largo = get<0>(camino[j]);
				int ancho = get<1>(camino[j]);
				if(_mapa_explorado[largo][ancho] == 0){

					_total_explorado +=1;
					if(objetivoCumplido()){
						return true;
					}
				}
				_mapa_explorado[largo][ancho] += 1;
			}
		}
		return false;
	}

	tupleList Solution::construirCamino(int vehiculo, bool greedy){
// printf("%s\n", "construirCamino");

		_movimientos_restantes_vehiculo[vehiculo] = _pbm.getAutonomiaVehiculo(vehiculo);
		// printf("%s", "_movimientos_restantes_vehiculo[vehiculo] ");
		// printf("%d\n", _movimientos_restantes_vehiculo[vehiculo]);
		// printf("%d\n", vehiculo);
		tupleList camino;
		// camino.emplace_back(0,0);
		bool finCamino = false;
		int largo = 0;
		int ancho = 0;
		_mapa_explorado[0][0] ++;
		bool retorno = false;
		bool first = true;
		while(! finCamino ){
			if(! retorno && !objetivoCumplido()){
				if (greedy){
					int minimo = 999999999;
					int mejorLargo = largo + 1;
					int mejorAncho = ancho + 1;
					for(int i = -1; i<2; i++){
						for(int j = -1; j <2; j++){
							int nuevoLargo = largo + i;
							int nuevoAncho = ancho + j;
							if ((!((i == j) && (i == 0))) && _pbm.isExplorable(nuevoLargo, nuevoAncho)){
								if ( (_mapa_explorado[nuevoLargo][nuevoAncho] < minimo)){
									minimo = _mapa_explorado[nuevoLargo][nuevoAncho];
									mejorLargo = nuevoLargo;
									mejorAncho = nuevoAncho;
								}
							}


						}
					}
					largo = mejorLargo;
					ancho = mejorAncho;
					explorarZona(largo, ancho, vehiculo);
					first = false;
						// camino.emplace_back(largo, ancho);
				}else{
					bool encontreLugar = false;
					while(! encontreLugar){
						int nuevoLargo = (rand() % 3) - 1;
						int nuevoAncho = (rand() % 3) - 1;
						nuevoLargo = largo + nuevoLargo;
						nuevoAncho = ancho + nuevoAncho;
						if (_pbm.isExplorable(nuevoLargo, nuevoAncho)){
							encontreLugar = true;
							largo = nuevoLargo;
							ancho = nuevoAncho;
							first = false;
							explorarZona(largo, ancho, vehiculo);
							// camino.emplace_back(largo, ancho);
						}
					}

				}





			}else{ //estoy volviendo a la base , FALTA ESQUIVAR OBSTACULOS

				if (largo != 0){

					largo = largo -1;
				}
				if (ancho != 0){
					ancho = ancho -1;
				}
				explorarZona(largo, ancho, vehiculo);

			}
			if (!first && _pbm.isBasePosition(largo, ancho)){
				// printf("%s\n", "isBasePosition");
				finCamino = true;
			}
			camino.emplace_back(largo, ancho);
			retorno = deboVolver(largo, ancho, vehiculo);


		}
		return camino;

	}

	bool Solution::objetivoCumplido() const
	{
		if(_total_explorado == 0 ){
			// printf("%s\n", "es 0");
			return false;
		}else{
			float porcentaje = ((double)_total_explorado / _pbm.getCantidadZonas());
			// printf("%s", "porcentaje ");
			// printf("%f\n", porcentaje);
			if (porcentaje >= 0.8){
				return true;
			} else{
				return false;
			}
		}
	}

	bool Solution::deboVolver(int largo, int ancho, int vehiculo){
		//si la cantidad de movimientos que debo realizar para volver a la base es igual a mi autonomia vuelvo
		bool finCamino = false;
		int distancia = largo;
		if (distancia > ancho){
			distancia = ancho;
		}
		if (distancia >= _movimientos_restantes_vehiculo[vehiculo]+1){
			finCamino =true;
		}
		return finCamino;
	}


	int Solution::explorarZona(int largo,int ancho, int vehiculo)
	{
		// printf("%d\n", largo);
		// printf("%d\n", ancho);
		// printf("%d\n", _mapa_explorado[largo][ancho]);
		if (_mapa_explorado[largo][ancho] == 0){
			_total_explorado ++;
			// printf("%s\n", "_total_explorado");
		}
		_mapa_explorado[largo][ancho] +=1;
		_movimientos_restantes_vehiculo[vehiculo] -=1;
		if(_pbm.isLoadZone(largo, ancho)){
			_movimientos_restantes_vehiculo[vehiculo] = _pbm.getAutonomiaVehiculo(vehiculo);
		}
		return _mapa_explorado[largo][ancho];
	}

	bool Solution::isAlreadyExplored(int largo,int ancho) const
	{
		return (_mapa_explorado[largo][ancho] > 0);
	}

	double Solution::fitness()
	{
		//recorremos cada camino y nos quedamos con el mas _largo
		int iterVehiculos = 0;
		int* largosCaminos = new int[_pbm.getCantidadVehiculos()];
		for(int j=0; j < _pbm.getCantidadVehiculos(); j++){
			largosCaminos[j] = 0;
		}
		for (int i=0;i<_var.size();i++){//recorro tupla solucion(tareas)
			if(_var[i] == -2){
				break;
			}
			if(_var[i] == -1){
				iterVehiculos ++;
				if (iterVehiculos == _pbm.getCantidadVehiculos()){
					iterVehiculos = 0;
				}
			}else{

				largosCaminos[iterVehiculos] += 1;
				i++;
			}

		}

		int largoMasLargo = largosCaminos[0];
		for(int j=1; j<_pbm.getCantidadVehiculos();j++){
			if(largosCaminos[j] > largoMasLargo){
				largoMasLargo = largosCaminos[j];
			}
		}

		return largoMasLargo;
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
		// Se implementa cruzamiento de un punto, que implica intercambiar las rutas para un vehiculo
		// Se controla que se mantenga una solucion factible
// printf("%s\n", "cross");
		sol1.convertVarToCaminos();
		sol2.convertVarToCaminos();
		vector<tupleList> caminosSol1 = sol1.getCaminos();
		vector<tupleList> caminosSol2 = sol2.getCaminos();
		// printf("%s\n", "CAMINO 1");
		// for (int i = 0; i < caminosSol1.size(); i ++){
		// 	tupleList camino = caminosSol1[i];
		// 	for(int j = 0; j < camino.size(); j++){
		// 		printf("%d-", get<0>(camino[j]));
		// 		printf("%d\n", get<1>(camino[j]));
		// 	}
		// 	printf("%s\n", "finCamino");
		// }
		// printf("%s\n", "CAMINO 2");
		// for (int i = 0; i < caminosSol2.size(); i ++){
		// 	tupleList camino = caminosSol2[i];
		// 	for(int j = 0; j < camino.size(); j++){
		// 		printf("%d-", get<0>(camino[j]));
		// 		printf("%d\n", get<1>(camino[j]));
		// 	}
		// 	printf("%s\n", "finCamino");
		// }
		// printf("%s\n", "camino 1");
		// for(int i = 0 ;i <sol1.array_var().size(); i++){
		// 	if(sol1.var(i) == -2){
		// 		break;
		// 	}
		// 		printf("%d-", sol1.var(i));
		// }
		// printf("%s\n", "camino 2");
		// for(int i = 0 ;i <sol2.array_var().size(); i++){
		// 	if(sol2.var(i) == -2){
		// 		break;
		// 	}
		// 		printf("%d-", sol2.var(i));
		// }
		int puntoDeCorte = rand_int(0,caminosSol1.size() -1);
		// printf("%s-", "puntoDeCorte");
		// printf("%d\n", puntoDeCorte);
		// for(int i = puntoDeCorte; i<caminosSol1.size(); i+=punto){
			if ((puntoDeCorte < caminosSol1.size())&& (puntoDeCorte < caminosSol2.size())){

				caminosSol1[puntoDeCorte].swap(caminosSol2[puntoDeCorte]);
			}

		// }
// printf("%s\n", "pasa for");
	// printf("%s", "porcentaje: ");
	// printf("%d\n", sol1.total_explorado());
	// printf("%d\n", sol1.pbm().getCantidadZonas());

		sol1.setCaminos(caminosSol1);
		sol2.setCaminos(caminosSol2);
		std::vector<tupleList>().swap(caminosSol1);
		std::vector<tupleList>().swap(caminosSol2);
		// printf("%s\n", "CAMINO 1");
		// for (int i = 0; i < caminosSol1.size(); i ++){
		// 	tupleList camino = caminosSol1[i];
		// 	for(int j = 0; j < camino.size(); j++){
		// 		printf("%d-", get<0>(camino[j]));
		// 		printf("%d\n", get<1>(camino[j]));
		// 	}
		// 	printf("%s\n", "finCamino");
		// }
		// printf("%s\n", "CAMINO 2");
		// for (int i = 0; i < caminosSol2.size(); i ++){
		// 	tupleList camino = caminosSol2[i];
		// 	for(int j = 0; j < camino.size(); j++){
		// 		printf("%d-", get<0>(camino[j]));
		// 		printf("%d\n", get<1>(camino[j]));
		// 	}
		// 	printf("%s\n", "finCamino");
		// }
		if(sol1.cumpleObjetivo() && sol2.cumpleObjetivo()){
			// printf("%d\n", ((double)sol1.total_explorado() / sol1.pbm().getCantidadZonas()));
			sol1.convertCaminosToVar();
			sol2.convertCaminosToVar();

			// printf("%s\n", "camino 1");
			// for(int i = 0 ;i <sol1.array_var().size(); i++){
			// 	if(sol1.var(i) == -2){
			// 		break;
			// 	}
			// 		printf("%d-", sol1.var(i));
			// }
			// printf("%s\n", "camino 2");
			// for(int i = 0 ;i <sol2.array_var().size(); i++){
			// 	if(sol2.var(i) == -2){
			// 		break;
			// 	}
			// 		printf("%d-", sol2.var(i));
			// }
		}
// printf("%s\n", "end cross");
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
		// printf("%s\n", "mutate");
		sol.convertVarToCaminos();
		// for(int i = 0 ;i <sol.array_var().size(); i++){
		// 	if(sol.var(i) == -2){
		// 		break;
		// 	}
		// 		printf("%d-", sol.var(i));
		// }
		vector<tupleList> caminosSol1 = sol.getCaminos();
		if(caminosSol1.size() >0){
			int caminoABorrar = rand_int(0,caminosSol1.size() -1);
			caminosSol1.erase(caminosSol1.begin() + caminoABorrar);
		// 	tupleList caminosVehiculo = caminosSol1[vehiculo];
		// // printf("%d\n", caminosSol1[vehiculo].size());
		// 	int caminoABorrar = rand_int(0,caminosVehiculo.size());
		// 	// printf("%d\n", caminosVehiculo.size());
		// 	printf("%s\n", "caminosVehicule");
		// 	for(int i = 0 ; i < caminosVehiculo.size(); i ++){
		// 		printf("%d-", get<0>(caminosVehiculo[i]));
		// 		printf("%d\n", get<1>(caminosVehiculo[i]));
		// 	}
		// 	caminosVehiculo.erase(caminosVehiculo.begin() + caminoABorrar);
		// 	 caminosSol1[vehiculo]= caminosVehiculo;
		// 	 printf("%s\n", "caminosVehiculeAfterremove");
		// 	 for(int i = 0 ; i < caminosVehiculo.size(); i ++){
		// 		 printf("%d-", get<0>(caminosVehiculo[i]));
		// 		 printf("%d\n", get<1>(caminosVehiculo[i]));
			sol.setCaminos(caminosSol1);

			std::vector<tupleList>().swap(caminosSol1);

			if(sol.cumpleObjetivo()){
				// printf("%s\n", "cumpleObjetivo");
				// printf("%s\n", "After mutate");
				sol.convertCaminosToVar();
				// for(int i = 0 ;i <sol.array_var().size(); i++){
				// 	if(sol.var(i) == -2){
				// 		break;
				// 	}
				// printf("%d-", sol.var(i));

			}
		}

			// printf("%d\n", caminosVehiculo.size());
			// printf("%d\n", caminosSol1[vehiculo].size());


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
