/* g++ makerandompath.cpp -o run
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iomanip> // precision //
#include <string.h>
#include <math.h>
#include <cstdlib> // random numbers //
#include <time.h> // for 'random' seed //

const int N = 512;
const int DIM = 2;
const double PI = 3.1415926535f;

class Point {
  public:
    double x, y;
    Point() {
      x = 0.0f;
      y = 0.0f;
    }
    Point(double input[DIM]) {
      x = input[0];
      y = input[1];
    }
    Point(double x0, double y0) {
      x = x0;
      y = y0;
    }
	 std::string str() {
		 	std::ostringstream tmp;
			tmp << std::string( "(" ) << x << std::string( ", " ) << y << std::string( ")" );
			return tmp.str();
	 }
};
class Path {
  public:
	 int npoints; // how many points are sampled?
	 int nmodes; // how many modes are taken into account? from [pi to nmodes*pi]
	 double amplitude; // amplitude of coefficients in path generation
	 double amplitude_limit; // limitation for values of rescaled amplitudes
	 int pathdensity; // number of samples per arc length 1
	 bool allocated; // test if array of points was allocated
    Point* points; // array of points
    Path() {
		nmodes = 5;
		amplitude = 1.0f;
		amplitude_limit = 1.0f;
		pathdensity = 100;
		allocated = false;
    }
    ~Path() {
      if (allocated) delete[] points;
    }
    void makepath_fourier () {
		// find cosine coefficients for x-coordinate
  		double x_coeffs_cos[nmodes], x_coeffs_sin[nmodes];
		bool amplitude_limit_satisfied = false;
		while (!amplitude_limit_satisfied) {
			// total sum = 0
			double thissum = 0.0f;
			for (int i = 2; i <= nmodes; i++) {
				x_coeffs_cos[i-1] = amplitude * 2.0f * (((double) (rand() % 10000)) / 10000.0f - 0.5f)/((double) (i*i));
				thissum += x_coeffs_cos[i-1];
			}
			x_coeffs_cos[0] = -1.0f;
	  		for (int i = 2; i <= nmodes; i++) {
	  			x_coeffs_cos[i-1] /= thissum;
	  		}
			// alternating sum = 1
			thissum = 0.0f;
			for (int i = 1; i <= nmodes; i++) {
				if (i % 2 == 0)
					thissum += x_coeffs_cos[i-1];
				else
					thissum -= x_coeffs_cos[i-1];
			}
			for (int i = 1; i <= nmodes; i++) {
				x_coeffs_cos[i-1] /= thissum;
			}
			// find sine coefficients
			for (int i = 1; i <= nmodes; i++) {
				x_coeffs_sin[i-1] = amplitude * 2.0f * (((double) (rand() % 10000)) / 10000.0f - 0.5f)/((double) (i*i));
			}
			amplitude_limit_satisfied = true;
			for (int i = 0; i < nmodes; i++) {
				if (fabs(x_coeffs_sin[i]) > amplitude_limit || fabs(x_coeffs_cos[i]) > amplitude_limit) {
					amplitude_limit_satisfied = false;
					break;
				}
			}
		}

  		// find cosine coefficients for y-coordinate
		double y_coeffs_cos[nmodes], y_coeffs_sin[nmodes];
		amplitude_limit_satisfied = false;
		while (!amplitude_limit_satisfied) {
  			// total sum = 0
			double thissum = 0.0f;
	  		for (int i = 2; i <= nmodes; i++) {
	  			y_coeffs_cos[i-1] = amplitude * 2.0f * (((double) (rand() % 10000)) / 10000.0f - 0.5f)/((double) (i*i));
	  			thissum += y_coeffs_cos[i-1];
	  		}
	  		y_coeffs_cos[0] = -1.0f;
	  		for (int i = 2; i <= nmodes; i++) {
	  			y_coeffs_cos[i-1] /= thissum;
	  		}
	  		// alternating sum = 1
	  		thissum = 0.0f;
	  		for (int i = 1; i <= nmodes; i++) {
	  			if (i % 2 == 0)
	  				thissum += y_coeffs_cos[i-1];
	  			else
	  				thissum -= y_coeffs_cos[i-1];
	  		}
	  		for (int i = 1; i <= nmodes; i++) {
	  			y_coeffs_cos[i-1] /= thissum;
	  		}
	  		// find sine coefficients
	  		for (int i = 1; i <= nmodes; i++) {
	  			y_coeffs_sin[i-1] = amplitude * 2.0f * (((double) (rand() % 10000)) / 10000.0f - 0.5f)/((double) (i*i));
	  		}
			amplitude_limit_satisfied = true;
			for (int i = 0; i < nmodes; i++) {
				if (fabs(y_coeffs_sin[i]) > amplitude_limit || fabs(y_coeffs_cos[i]) > amplitude_limit) {
					amplitude_limit_satisfied = false;
					break;
				}
			}
		}

/*
		// output amplitudes
  		for (int i = 0; i < nmodes; i++) {
			std::cout << "x at " << i+1 << "pi: " << x_coeffs_cos[i] << " " << x_coeffs_sin[i] << std::endl;
  		}
		std::cout << "=============" << std::endl;
	  	for (int i = 0; i < nmodes; i++) {
			std::cout << "y at " << i+1 << "pi: " << y_coeffs_cos[i] << " " << y_coeffs_sin[i] << std::endl;
  		}
*/

		// estimate length of path
		double path_arclength = 0.0f;
		double currentt = 0.0f, thisx = 0.0f, thisy = 0.0f, nextx = 0.0f, nexty = 0.0f;
  		for (int i = 0; i < 200; i++) {
			nextx = thisx;
			nexty = thisy;
			nextx = calcpath(x_coeffs_cos, x_coeffs_sin, currentt);
			nexty = calcpath(y_coeffs_cos, y_coeffs_sin, currentt);
			path_arclength += sqrt((thisx - nextx)*(thisx - nextx) + (thisy - nexty)*(thisy - nexty));
			currentt = (double) (i + 1) / 200.0f;
  		}
		allocatepath(floor(path_arclength * pathdensity + 0.5f));

		for (int i = 0; i < npoints; i++) {
			currentt = (double) i / (double) (npoints - 1);
			points[i].x = calcpath(x_coeffs_cos, x_coeffs_sin, currentt);
			points[i].y = calcpath(y_coeffs_cos, y_coeffs_sin, currentt);
		}
    }
  private:
	  double calcpath (double* coeffs_cos, double* coeffs_sin, double t) {
		  double result = 0.0f;
		  for (int i = 1; i <= nmodes; i++) {
			  result += coeffs_cos[i - 1] * cos(i * PI * t) + coeffs_sin[i - 1] * sin(i * PI * t);
		  }
		  return result;
	  }
	  void allocatepath (int npoints0) {
			if (allocated) delete[] points;
			npoints = npoints0;
			points = new Point[npoints];
			allocated = true;
	  }
};

void writepathtofile( Path somepath );

int main(int argc, char *argv[]){
	unsigned int someseed = (unsigned int) (time(NULL));
	srand(someseed);

	Path somepath;
	somepath.makepath_fourier();
	writepathtofile(somepath);

	return 0;
}

void writepathtofile( Path somepath ){
	std::ofstream outputfile;
	int icount;

	outputfile.open("output.dat");
	for (icount = 0; icount < somepath.npoints; icount++){
		outputfile << somepath.points[icount].x << "\t" << somepath.points[icount].y << std::endl;
	}
	outputfile.close();

	return;
}
