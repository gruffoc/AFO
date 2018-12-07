#include <vector>

class Integrate {

	private:
		double (*fun_)(double, double, double);
		std::vector< double > _ranges;
		int _dim;
		int _MCSteps;
		


	public:
		Integrate(double (*fun)(double, double, double), std::vector< double > ranges, int MCSteps);  //CTOR
		~Integrate(); //DTOR
		
		double cubatura_naif();


};
