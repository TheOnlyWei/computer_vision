// By Wei Shi (based on code by Ioannis Stamos)
// last modified 11/10/16
// Class for representing a 2D gray-scale image,
// with support for reading/writing pgm images.


#include "image.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

namespace ComputerVisionProjects {

Image::Image(const Image &an_image){

  AllocateSpaceAndSetSize(an_image.num_rows(), an_image.num_columns());
  SetNumberGrayLevels(an_image.num_gray_levels());
	labels = an_image.labels;
	//need_update = an_image.need_update;


  for (size_t i = 0; i < num_rows(); ++i)
    for (size_t j = 0; j < num_columns(); ++j){
      SetPixel(i,j, an_image.GetPixel(i,j));
    }
}

Image& Image::operator=(const Image &an_image) {

  AllocateSpaceAndSetSize(an_image.num_rows(), an_image.num_columns());
  SetNumberGrayLevels(an_image.num_gray_levels());
	labels = an_image.labels;
	//need_update = an_image.need_update;

  for (size_t i = 0; i < num_rows(); ++i)
    for (size_t j = 0; j < num_columns(); ++j){
      SetPixel(i,j, an_image.GetPixel(i,j));
    }

	return *this;

}

Image::~Image(){
  DeallocateSpace();
}

void Image::AllocateSpaceAndSetSize(size_t num_rows, size_t num_columns) {

  if (pixels_ != nullptr) {
		DeallocateSpace();
	}

  pixels_ = new int*[num_rows];

  for (size_t i = 0; i < num_rows; ++i)
    pixels_[i] = new int[num_columns];

  num_rows_ = num_rows;
  num_columns_ = num_columns;
}
void Image::allocateFlatGrayScale(size_t num_rows, size_t num_columns, unsigned short gray_level) {

  if (pixels_ != nullptr) {
		DeallocateSpace();
	}

  pixels_ = new int*[num_rows];

  for (size_t i = 0; i < num_rows; ++i){
    pixels_[i] = new int[num_columns];

		for(size_t j = 0; j < num_columns; ++j) {

			pixels_[i][j] = gray_level;
		}
	}

  num_rows_ = num_rows;
  num_columns_ = num_columns;


}
void Image::DeallocateSpace() {

  for (size_t i = 0; i < num_rows_; i++) {
    delete pixels_[i];
	}
  delete pixels_;

  pixels_ = nullptr;
  num_rows_ = 0;
  num_columns_ = 0;
}


void Image::printGrayscales() {

	cout << setw(5);
	for(size_t i = 0; i < num_rows(); ++i) {
		for(size_t j = 0; j < num_columns(); ++j) {
			if(i > 0 && j % num_columns() == 0)
				cout << endl;

			cout << pixels_[i][j] << setw(5);

		}
	}

	cout << endl;
}


/*
unordered_map<unsigned short, double> Image::getCenterX() {

	if(need_update) {
		setCenterX();
	}
	return centerX;
}

unordered_map<unsigned short, double> Image::getCenterY() {
	if(need_update) {
		setCenterY();
	}
	return centerY;

}


void Image::update() {
	need_update = true;

}
*/
/////////////////////////////////////////////////////////////////////////////////////////////
//
//																NON MEMBER FUNCTIONS
//
/////////////////////////////////////////////////////////////////////////////////////////////



unordered_map<unsigned short, double> getCenterX(Image* image) {

	/*
		For each row, iterate through the columns
		multiply row number by the binary pixel of the object
		sum each of such pixels

	*/

	unordered_map<unsigned short, unsigned long> label_to_area = getArea(image);
	unordered_map<unsigned short, double> label_to_x_moment;
	unordered_map<unsigned short, double>::const_iterator it;

	int pixel;

	for(size_t i = 0; i < image->num_rows(); ++i) {

		for(size_t j = 0; j < image->num_columns(); ++j) {

			pixel = image->GetPixel(i, j);
			it = label_to_x_moment.find(pixel);

			if(pixel != 0) {

				if(it == label_to_x_moment.end() ) {
					label_to_x_moment[pixel] = i;
				}
				else {
					label_to_x_moment[pixel] += i;

				}

			}
		}
	}

	unsigned long area;
	for(auto& pair: label_to_x_moment) {
		area = label_to_area[pair.first];
		pair.second /= area;

	}

	return label_to_x_moment;
}

unordered_map<unsigned short, double> getCenterY(Image* image) {

	unordered_map<unsigned short, unsigned long> label_to_area = getArea(image);
	unordered_map<unsigned short, double> label_to_y_moment;
	unordered_map<unsigned short, double>::const_iterator it;

	int pixel;

	for(size_t i = 0; i < image->num_rows(); ++i) {

		for(size_t j = 0; j < image->num_columns(); ++j) {

			pixel = image->GetPixel(i, j);
			it = label_to_y_moment.find(pixel);

			if(pixel != 0) {

				if(it == label_to_y_moment.end() ) {
					label_to_y_moment[pixel] = j;
				}
				else {
					label_to_y_moment[pixel] += j;

				}

			}


		}
	}

	unsigned long area;
	for(auto& pair: label_to_y_moment) {
		area = label_to_area[pair.first];
		pair.second /= area;

	}

	return label_to_y_moment;


}
bool ReadImage(const string &filename, Image *an_image) {
  if (an_image == nullptr) abort();

	//an_image->update();
	// fopen opens file as text files by default
	// to have fopen open files as binary, append a "b" to the second parameter (mode)
  FILE *input = fopen(filename.c_str(),"rb");

	// nullptr has value 0
  if (input == 0) {
    cout << "ReadImage: Cannot open file" << endl;
    return false;
  }

  // Check for the right "magic number".
  char line[1024];
	// size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );
	// Reads an array of count elements, each one with a size of size bytes,
	// from the stream and stores them in the block of memory specified by ptr.
  if (fread(line, 1, 3, input) != 3 || strncmp(line,"P5\n",3)) {
    fclose(input);
    cout << "ReadImage: Expected .pgm file" << endl;
    return false;
  }

  // Skip comments.
  do {
		// char * fgets ( char * str, int num, FILE * stream );
		// Reads characters from stream and stores them as a C string into str
		// until (num-1) characters have been read or either a newline
		// or the end-of-file is reached, whichever happens first.
		// sizeof returns size in bytes of an object
    fgets(line, sizeof line, input);

	}

  while(*line == '#');

  // Read the width and height.
  int num_columns,num_rows;
  sscanf(line,"%d %d\n", &num_columns, &num_rows);
  an_image->AllocateSpaceAndSetSize(num_rows, num_columns);


  // Read # of gray levels.
  fgets(line, sizeof line, input);
  int levels;
  sscanf(line,"%d\n", &levels);
  an_image->SetNumberGrayLevels(levels);

  // read pixel row by row.
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0;j < num_columns; ++j) {
      const int byte=fgetc(input);
      if (byte == EOF) {
        fclose(input);
        cout << "ReadImage: short file" << endl;
        return false;
      }
      an_image->SetPixel(i, j, byte);
    }
  }

  fclose(input);
  return true;
}

bool WriteImage(const string &filename, const Image &an_image) {
  FILE *output = fopen(filename.c_str(), "w");
  if (output == 0) {
    cout << "WriteImage: cannot open file" << endl;
    return false;
  }
  const int num_rows = an_image.num_rows();
  const int num_columns = an_image.num_columns();
  const int colors = an_image.num_gray_levels();

  // Write the header.
  fprintf(output, "P5\n"); // Magic number.
  fprintf(output, "#\n");  // Empty comment.
  fprintf(output, "%d %d\n%03d\n", num_columns, num_rows, colors);

  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_columns; ++j) {
      const int byte = an_image.GetPixel(i , j);
      if (fputc(byte,output) == EOF) {
	    fclose(output);
            cout << "WriteImage: could not write" << endl;
	    return false;
      }
    }
  }

  fclose(output);
  return true;
}

// Implements the Bresenham's incremental midpoint algorithm;
// (adapted from J.D.Foley, A. van Dam, S.K.Feiner, J.F.Hughes
// "Computer Graphics. Principles and practice",
// 2nd ed., 1990, section 3.2.2);
void DrawLine(int x0, int y0, int x1, int y1, int color, Image *an_image) {
  if (an_image == nullptr) abort();

#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) {a^=b; b^=a; a^=b;}

  const int DIR_X = 0;
  const int DIR_Y = 1;

  // Increments: East, North-East, South, South-East, North.
  int incrE,
    incrNE,
    incrS,
    incrSE,
    incrN;
  int d;         /* the D */
  int x,y;       /* running coordinates */
  int mpCase;    /* midpoint algorithm's case */
  int done;      /* set to 1 when done */

  int xmin = x0;
  int xmax = x1;
  int ymin = y0;
  int ymax = y1;

  int dx = xmax - xmin;
  int dy = ymax - ymin;
  int dir;

  if (dx * dx > dy * dy) {  // Horizontal scan.
    dir=DIR_X;
    if (xmax < xmin) {
      SWAP(xmin, xmax);
      SWAP(ymin , ymax);
    }
    dx = xmax - xmin;
    dy = ymax - ymin;

    if (dy >= 0) {
      mpCase = 1;
      d = 2 * dy - dx;
    } else {
      mpCase = 2;
      d = 2 * dy + dx;
    }

    incrNE = 2 * (dy - dx);
    incrE = 2 * dy;
    incrSE = 2 * (dy + dx);
  } else {// vertical scan.
    dir = DIR_Y;
    if (ymax < ymin) {
      SWAP(xmin, xmax);
      SWAP(ymin, ymax);
    }
    dx = xmax - xmin;
    dy = ymax-ymin;

    if (dx >=0 ) {
      mpCase = 1;
      d = 2 * dx - dy;
    } else {
      mpCase = 2;
      d = 2 * dx + dy;
    }

    incrNE = 2 * (dx - dy);
    incrE = 2 * dx;
    incrSE = 2 * (dx + dy);
  }

  /// Start the scan.
  x = xmin;
  y = ymin;
  done = 0;

  while (!done) {
    an_image->SetPixel(x,y,color);

    // Move to the next point.
    switch(dir) {
    	case DIR_X:
      	if (x < xmax) {
			    switch(mpCase) {
					  case 1:
						if (d <= 0) {
							d += incrE;
							x++;
						} else {
							d += incrNE;
							x++;
							y++;
						}
						break;

            case 2:
              if (d <= 0) {
                d += incrSE;
								x++;
								y--;
              } else {
                d += incrE;
								x++;
              }
					  break;
					 }
      } else {
				done=1;
      }
      break;

    case DIR_Y:
      if (y < ymax) {
        switch(mpCase) {
					case 1:
						if (d <= 0) {
							d += incrE;
							y++;
						} else {
							d += incrNE;
							y++;
							x++;
						}
					 break;

					case 2:
						if (d <= 0) {
							d += incrSE;
							y++;
							x--;
						} else {
						  d += incrE;
							y++;
						}
						break;
	  		} // mpCase
			} // y < ymax
      else {
	  		done=1;
			}
			break;
    } // switch(DIR)
  } // while(!done)
}

void toBinary(Image* output, unsigned short threshold) {
	// for each pixel
	// make it white
	// if p[i] < threshold => p[i] = 0
	// else p[i] = 1 (when p[i] >= threshold)

	for(size_t i = 0; i < output->num_rows(); ++i) {

		for(size_t j = 0; j < output->num_columns(); ++j) {
			if(output->GetPixel(i, j) < threshold) {
				output->SetPixel(i, j, 0);
			} else {
				output->SetPixel(i, j, 255);

			}
		}

	}

}

DisjointSets label(Image* output, unsigned short min, unsigned short max, size_t label_size) {

	int left = -1, up_left = -1, up = -1;
	unsigned short set_size;

	DisjointSets labels(label_size);
	vector<bool> valid_labels(label_size, false);
	unsigned short count = 1;

	// initialize p[0][0]
	if(output->GetPixel(0, 0) == max) {
		output->SetPixel(0, 0, count);
		valid_labels[count] = true;
		count++;
		if(count > label_size) {
			throw "ERROR: Number of labels are too small!";

		}
	}

	// initialize first row after column 0
	for(size_t j = 1; j < output->num_columns(); ++j) {
		int pre_pixel_shade = output->GetPixel(0, j-1);
		int cur_pixel_shade = output->GetPixel(0, j);

		if(cur_pixel_shade == max) {
			if(pre_pixel_shade > min) { // > 0 means not black or is valid
				output->SetPixel(0, j, pre_pixel_shade);
			}
			else {
				output->SetPixel(0, j, count);
				valid_labels[count] = true;
				count++;
				if(count > label_size) {
					throw "ERROR: Number of labels are too small!";

				}
			}

		}
	}



	// initialize first column after row 0
	for(size_t i = 1; i < output->num_rows(); ++i) {
		int pre_pixel_shade = output->GetPixel(i-1, 0);
		int cur_pixel_shade = output->GetPixel(i, 0);

		if(cur_pixel_shade == max) {
			if(pre_pixel_shade > min) { // > 0 means not black or is valid
				output->SetPixel(i, 0, pre_pixel_shade);
			}
			else {
				output->SetPixel(i, 0, count);
				valid_labels[count] = true;
				count++;
				if(count > label_size) {
					throw "ERROR: Number of labels are too small!";

				}
			}

		}
	}



	// loop through the rest
	for(size_t i = 1; i < output->num_rows(); ++i) {

		//cout << endl;
		for(size_t j = 1; j < output->num_columns(); ++j) {
			if(output->GetPixel(i, j) == max) {
				// find neighboring shade
				left = output->GetPixel(i, j-1);
				up_left = output->GetPixel(i-1, j-1);
				up = output->GetPixel(i-1, j);

				// assign new shade since no former neighbors are detected
				if(left <= min && up_left <= min && up <= min)  {
					output->SetPixel(i, j, count);
					valid_labels[count] = true;
					count++;

					if(count > label_size) {
						throw "ERROR: Number of labels are too small!";

					}

				}

				// two pixels match
				// left and up_left match
				else if(left > min && left == up_left) {
					output->SetPixel(i, j, up_left);
				}
				// up_left and up match
				else if(up_left > min && up_left == up) {
					output->SetPixel(i, j, up);
				}
				// up and left match
				else if(up > min && up == left) {
					output->SetPixel(i, j, up);
				}
				// previous pixel shades all match
				else if(left > min && left == up_left && up_left == up) {
					output->SetPixel(i, j, up);
				}

				// up and left are joined by pixel[i][j], but have different values
				else if (up > min && left > min && up != left) {
					// invariant: up, left, and up_left will always have their labels
					labels.unionSets( labels.find(up), labels.find(left) );

					output->SetPixel(i, j, up);
				}
				else if (up > min && up_left > min && up != up_left) {
					labels.unionSets( labels.find(up), labels.find(up_left) );

					output->SetPixel(i, j, up_left);
				}
				else if (left > min && up_left > min && left != up_left) {
					labels.unionSets( labels.find(left), labels.find(up_left) );

					output->SetPixel(i, j, up_left);
				}

				else if(up_left > min) {
					output->SetPixel(i, j, up_left);

				}

				else if(up > min) {
					output->SetPixel(i, j, up);
				}

				else if(left > min) {
					output->SetPixel(i, j, left);
				}

			}
		}

	}

	labels.setValidIndex(valid_labels);

	return labels;
}



unordered_map<size_t, unsigned short> generateShadeMap(const DisjointSets& labels,
																													 unsigned short base,
																													 short increment) {

	unordered_map<size_t, unsigned short> label_to_gray;

	const vector<int> label_set = labels.getDisjointSet();

	// set all roots first
	for(size_t i = 0; i < label_set.size(); ++i) {

		unordered_map<size_t, unsigned short>::const_iterator it;
		it = label_to_gray.find(labels.find(i));

		if( it == label_to_gray.end() && labels.isUsed(i) ) {
			// initialize root's shade

			int root = labels.find(i);
			label_to_gray[root] = base;
			base += increment;
		}

	}

	for(size_t i = 0; i < label_set.size(); ++i) {

		unordered_map<size_t, unsigned short>::const_iterator it;
		it = label_to_gray.find(i);

		if( it == label_to_gray.end() && labels.isUsed(i) ) {
			int root = labels.find(i);
			int root_shade = label_to_gray[root];
			label_to_gray[i] = label_to_gray[root];

		}

	}

	return label_to_gray;

}




void separateObjectsByShade(Image* image, unsigned short base,
																					short increment,
																					unsigned short min,
																					unsigned short max,
																					size_t label_size) {


	DisjointSets labels;
	labels = label(image, min, max, label_size);

	// for each i j in label_set and pixel array
	// we look to see if pixel[i][j] and pixel[i][j+1] have the same root
	// if they have the same root, we paint them the same grayscale
	unordered_map<size_t, unsigned short> label_to_gray = generateShadeMap(labels, base, increment);

	for(size_t i = 0; i < image->num_rows(); ++i) {

		for(size_t j = 0; j < image->num_columns(); ++j) {

			size_t pixel_value = image->GetPixel(i, j);

			if(pixel_value == 0)
				continue;

			unsigned short shade = label_to_gray[pixel_value];
			image->SetPixel(i, j, shade);

		}

	}

}


unordered_map<unsigned short, unsigned long> getArea(Image* image) {

	unordered_map<unsigned short, unsigned long> label_to_area;
	unordered_map<unsigned short, unsigned long>::const_iterator it;
	int pixel;

	for(size_t i = 0; i < image->num_rows(); ++i) {
		for(size_t j = 0; j < image->num_columns(); ++j) {

			pixel = image->GetPixel(i, j);

			if(pixel != 0) {
				it = label_to_area.find(pixel);

				if(it == label_to_area.end()) {
					label_to_area[pixel] = 1;
				}
				else {
					label_to_area[pixel]++;
				}

			}


		} // end column

	} // end row

	return label_to_area;

}

unordered_map<unsigned short, unsigned long> getNonhomogenousArea(Image* labeled_image,
																																	Image* nonhomogenous_image) {

	unordered_map<unsigned short, unsigned long> label_to_area;
	unordered_map<unsigned short, unsigned long>::const_iterator it;
	int pixel;

	for(size_t i = 0; i < labeled_image->num_rows(); ++i) {
		for(size_t j = 0; j < labeled_image->num_columns(); ++j) {

			pixel = labeled_image->GetPixel(i, j);

			if(pixel != 0) {
				it = label_to_area.find(pixel);

				if(it == label_to_area.end()) {
					label_to_area[pixel] = nonhomogenous_image->GetPixel(i, j);
				}
				else {
					label_to_area[pixel] += nonhomogenous_image->GetPixel(i, j);
				}

			}


		} // end column

	} // end row


	return label_to_area;

}


unordered_map<unsigned short, vector<double>> getSecondMoments(Image* input) {
	// v[0] = a, v[1] = b, v[2] = c
	unordered_map<unsigned short, vector<double> > label_to_2nd_moments;
	unordered_map<unsigned short, vector<double> >::const_iterator it;

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(input);
	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(input);

	int pixel;
	double x_center, y_center;
	double x_prime, y_prime;
	vector<double> cur_moments;

	for(size_t i = 0; i < input->num_rows(); ++i) {

		for(size_t j = 0; j < input->num_columns(); ++j) {
			pixel = input->GetPixel(i, j);

			if(pixel != 0) {
				it = label_to_2nd_moments.find(pixel);
				x_center = label_to_x_moment[pixel];
				y_center = label_to_y_moment[pixel];
				x_prime = i - x_center, y_prime = j - y_center;

				if(it == label_to_2nd_moments.end()) {
					cur_moments = vector<double>();
					cur_moments.push_back(x_prime*x_prime);
					cur_moments.push_back(x_prime*y_prime);
					cur_moments.push_back(y_prime*y_prime);
					label_to_2nd_moments[pixel] = cur_moments;

				} else {
					label_to_2nd_moments[pixel][0] += x_prime * x_prime;
					label_to_2nd_moments[pixel][1] += x_prime * y_prime;
					label_to_2nd_moments[pixel][2] += y_prime * y_prime;

				}
			}
		}

	}

	// Multiply moment 'b' by 2
	for(auto& pair: label_to_2nd_moments) {
		pair.second[1] *= 2;
	}

	return label_to_2nd_moments;
}

unordered_map<unsigned short, double> getTheta(Image* input) {

	unordered_map<unsigned short, vector<double> > label_to_2nd_moments = getSecondMoments(input);
	unordered_map<unsigned short, double> label_to_theta;

	double theta;
	unsigned short label;
	long rhs_denom;

	long a, b, c;

	for(const auto& pair: label_to_2nd_moments) {
		label = pair.first;
		a = pair.second[0], b = pair.second[1], c = pair.second[2];
		rhs_denom = a - c;
		theta = ( atan2(b, rhs_denom) ) / 2;
		label_to_theta[label] = theta;

	}

	return label_to_theta;
}

unordered_map<unsigned short, double> getEmin(Image* input) {

	unordered_map<unsigned short, double> label_to_theta = getTheta(input);
	unordered_map<unsigned short, vector<double>> label_to_2nd_moments = getSecondMoments(input);
	unordered_map<unsigned short, double> label_to_Emin;

	double sin_theta, cos_theta, Emin;
	unsigned short label;
	double a, b, c;

	for(const auto& pair: label_to_theta) {
		label = pair.first;

		a = label_to_2nd_moments[label][0];
		b = label_to_2nd_moments[label][1];
		c = label_to_2nd_moments[label][2];

		sin_theta = sin(pair.second);
		cos_theta = cos(pair.second);

		Emin = (a*pow(sin_theta, 2)) - (b * sin_theta * cos_theta) + (c * pow(cos_theta, 2));

		label_to_Emin[label] = Emin;

	}

	return label_to_Emin;
}

unordered_map<unsigned short, double> getEmax(Image* input) {

	unordered_map<unsigned short, double> label_to_theta = getTheta(input);
	unordered_map<unsigned short, vector<double>> label_to_2nd_moments = getSecondMoments(input);
	unordered_map<unsigned short, double> label_to_Emax;

	double sin_2theta, cos_2theta, Emax;
	double sin_theta, cos_theta, theta;
	unsigned short label;
	double a, b, c;


	for(const auto& pair: label_to_theta) {
		label = pair.first;
		theta = (M_PI/2) + pair.second;

		a = label_to_2nd_moments[label][0];
		b = label_to_2nd_moments[label][1];
		c = label_to_2nd_moments[label][2];

		sin_theta = sin(theta);
		cos_theta = cos(theta);

		Emax = (a*pow(sin_theta, 2)) - (b * sin_theta * cos_theta) + (c * pow(cos_theta, 2));

		label_to_Emax[label] = Emax;
	}

	return label_to_Emax;
}

unordered_map<unsigned short, double> getRoundness(Image* input) {

	unordered_map<unsigned short, double> label_to_Emin = getEmin(input);
	unordered_map<unsigned short, double> label_to_Emax = getEmax(input);
	unordered_map<unsigned short, double> label_to_roundness;

	double roundness, Emax, Emin;
	unsigned short label;

	for(const auto& pair: label_to_Emin) {
		label = pair.first;
		Emax = label_to_Emax[label];
		Emin = label_to_Emin[label];

		if(Emax != 0)
			roundness = label_to_Emin[label] / Emax;

		else if(Emin == Emax)
			roundness = 1;

		label_to_roundness[label] = roundness;
	}

	return label_to_roundness;
}

double gaussian(const double& mu, const double& sigma, const double& target) {
	double exp_numerator = -1 *pow( (target - mu), 2);
	double exp_denominator = 2 * pow(sigma, 2);
	double exp_ = exp(exp_numerator/ exp_denominator);
	double divisor = sqrt(2.0 * M_PI * pow(sigma, 2));

	return (exp_ / divisor);

}


void drawMoments(Image* input, Image* output, unsigned short dot_size) {

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(input);
	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(input);

	size_t half_dot;
	unsigned long x_moment;
	unsigned long y_moment;
	unsigned short label;

	for(const auto& pair: label_to_x_moment) {
		label = pair.first;
		//dot_size = label_to_dot_size[label];
		half_dot = dot_size/2;

		x_moment = round(pair.second);
		y_moment = round(label_to_y_moment.at(label));

		if(half_dot == 0) {

			output->SetPixel(x_moment, y_moment, 255);
		}

		else {
			for(size_t i = (x_moment - half_dot); i < (x_moment + half_dot) &&
																						i < output->num_rows(); ++i) {

				for(size_t j = (y_moment - half_dot); j < (y_moment + half_dot) &&
																							j < output->num_columns(); ++j) {

					output->SetPixel(i, j, 255);
				}

			}

		}

	}

}


void drawOrientation(Image* input, Image* output) {

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(input);
	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(input);
	unordered_map<unsigned short, vector<double> > label_to_2nd_moments = getSecondMoments(input);
	unordered_map<unsigned short, double> label_to_theta = getTheta(input);
	unordered_map<unsigned short, double> label_to_Emin = getEmin(input);

	// one point is at (v[0], v[1]), the other is at (v[2], v[3])
	unordered_map<unsigned short, vector<unsigned long> > label_to_points;

	double tan_theta, theta, constant;
	size_t i, j, last_x, last_y;
	unsigned short label;
	double x_center, y_center;
	size_t pixel;

	unsigned long radius;
	double distance;
	double lower_bound = (M_PI/2) - 0.0001, upper_bound = (M_PI/2) + 0.0001;

	vector<unsigned long> points(4, 0);

	/*
		The alternative algorithm is to use the formula:
		x * sin(theta) - y * cos(theta) + p = 0
		To find the equation, we need to solve for p,
		plug in x = x_bar, y = y_bar:

		x_bar * sin(theta) - y_bar * cos(theta) + p = 0

		Since theta can be calculated from the second moments, then
		the only unknown is p:

		-x_bar * sin(theta) + y_bar * cos(theta) = p
	*/


	// this drawing algorithm doesn't work well for sparse/noisy objects
	for(const auto& pair: label_to_theta) {

		label = pair.first;
		theta = pair.second;

		x_center = round(label_to_x_moment[label]);
		y_center = round(label_to_y_moment[label]);

		points[0] = x_center;
		points[1] = y_center;


		// x in terms of y
		if(lower_bound < theta && theta < upper_bound) {

			tan_theta = 0;
			constant = 	label_to_x_moment[label];
			i = x_center;
			last_x = x_center;
			last_y = y_center;
			j = y_center -1;

			// travel all the way to the edge of the image plane
			// record last seen pixel with current label.
			while(j >= 0 && j < input->num_columns()) {
				//distance = sqrt( (x_center - i)*(x_center - i) + (y_center - j)*(y_center - j) );
				// i doesn't change since it is a line parallel to y-axis
				pixel = input->GetPixel(i, j);

				// record last known label coordinate
				if(pixel == label) {
					// no need to check i since it stays the same.
					last_y = j;
				}

				if(j != 0)
					j--;
				else
					break;
			}

		}
		// y in terms of x
		else {
			tan_theta = tan(theta);
			constant = 	(-1 * label_to_x_moment[label] * tan_theta) + label_to_y_moment[label];
			//trace to one edge on axis of least inertia

			i = (x_center == 0 ) ? 0 : x_center - 1;
			j = y_center;
			last_x = x_center;
			last_y = y_center;

			while(i >= 0 && j>=0 && i < input->num_rows() && j < input->num_columns()) {

				if(round( (tan_theta * i) + constant )  < input->num_columns() &&
					 round( (tan_theta * i) + constant )  >= 0) {

					j = round( (tan_theta * i) + constant );
					pixel = input->GetPixel(i, j);
				}

				if(pixel == label) {
					last_x = i;
					last_y = j;

				}

				if(i != 0) {

					i--;
				}
				else
					break;
			}
		}
		points[2] = last_x;
		points[3] = last_y;

		label_to_points[label] = points;
		points = vector<unsigned long>(4,0);
	}

	int x0, y0, x1, y1;
	for(const auto& pair: label_to_points) {
		x0 = pair.second[0], y0 = pair.second[1], x1 = pair.second[2], y1 = pair.second[3];
		DrawLine(x0, y0, x1, y1, 255, output);
	}


}

void drawSelectLine(Image* input, Image* output, vector<unsigned short> labels) {

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(input);
	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(input);
	unordered_map<unsigned short, double> label_to_theta = getTheta(input);
	unordered_map<unsigned short, vector<unsigned long> > label_to_points;

	double lower_bound = (M_PI/2) - 0.0001, upper_bound = (M_PI/2) + 0.0001;
	double x_center, y_center;
	double tan_theta, theta, constant;
	size_t i, j, last_x, last_y;
	size_t pixel;

	vector<unsigned long> points(4, 0);

	for(const auto& label: labels) {

		x_center = round(label_to_x_moment[label]);
		y_center = round(label_to_y_moment[label]);
		theta = label_to_theta[label];
		points[0] = x_center;
		points[1] = y_center;

		if(lower_bound < theta && theta < upper_bound) {
			tan_theta = 0;
			constant = 	label_to_x_moment[label] ;

			i = x_center;
			last_x = x_center;
			last_y = y_center;
			j = y_center -1;

			// travel all the way to the edge of the image plane
			// record last seen pixel with current label.
			while(j >= 0 && j < input->num_columns()) {
				//distance = sqrt( (x_center - i)*(x_center - i) + (y_center - j)*(y_center - j) );
				// i doesn't change since it is a line parallel to y-axis
				pixel = input->GetPixel(i, j);

				// record last known label coordinate
				if(pixel == label) {
					// no need to check i since it stays the same.
					last_y = j;
				}

				if(j != 0)
					j--;
				else
					break;
			}

		}
		// y in terms of x
		else {

			tan_theta = tan(theta);
			constant = 	(-1 * label_to_x_moment[label] * tan_theta) + label_to_y_moment[label];
			//trace to one edge on axis of least inertia

			i = (x_center == 0 ) ? 0 : x_center - 1;
			j = y_center;
			last_x = x_center;
			last_y = y_center;

			while(i >= 0 && j>=0 && i < input->num_rows() && j < input->num_columns()) {

				if(round( (tan_theta * i) + constant )  < input->num_columns() &&
					 round( (tan_theta * i) + constant )  >= 0) {

					j = round( (tan_theta * i) + constant );
					pixel = input->GetPixel(i, j);
				}

				if(pixel == label) {
					last_x = i;
					last_y = j;

				}

				if(i != 0) {

					i--;
				}
				else
					break;
			}
		} // else
		points[2] = last_x;
		points[3] = last_y;

		label_to_points[label] = points;
		points = vector<unsigned long>(4,0);
	}

	int x0, y0, x1, y1;

	for(const auto& pair: label_to_points) {
		x0 = pair.second[0], y0 = pair.second[1], x1 = pair.second[2], y1 = pair.second[3];
		DrawLine(x0, y0, x1, y1, 255, output);
	}

}

void drawSelectCenter(Image* input, Image* output, vector<unsigned short>labels, size_t dot_size) {

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(input);
	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(input);

	size_t half_dot = dot_size/2;
	unsigned long x_center;
	unsigned long y_center;

	for(const auto& label: labels) {
		x_center = round(label_to_x_moment.at(label));
		y_center = round(label_to_y_moment.at(label));

		if(half_dot == 0) {
			output->SetPixel(x_center, y_center, 255);
		}

		else {
			for(size_t i = (x_center - half_dot); i < (x_center + half_dot) &&
																						i < output->num_rows(); ++i) {

				for(size_t j = (y_center - half_dot); j < (y_center + half_dot) &&
																							j < output->num_columns(); ++j) {

					output->SetPixel(i, j, 255);
				}

			}

		}

	} // end for

}

vector<vector<int>> deriveGaussian(const vector<vector<int>>& input, vector<unsigned int> kernel) {
	vector<vector<int>> temp_vector(2, vector<int>(2));
	vector<vector<int>> output(2, vector<int>(2));
	const size_t kernel_width = kernel.size();
	long product;
	long convolve_sum = 0;
	const int left = -1 * (kernel_width / 2);

	long temp;

	for(size_t i = 0; i < input.size(); ++i) {
		for(size_t j = 0; j < input[0].size(); ++j) {
			int k = left;
			for(size_t kernel_index = 0; kernel_index < kernel_width; ++kernel_index, ++k) {
				temp = j - k;
				if(temp < input[0].size() && temp >= 0 ) {
					product = (input[i][temp] * kernel[kernel_index]);
					convolve_sum += product;
				}
			}
			temp_vector[i][j] = convolve_sum;
			convolve_sum = 0, temp = 0;
		}

	}


	for(size_t j = 0; j < temp_vector.size(); ++j) {
		for(size_t i = 0; i < temp_vector[0].size(); ++i) {
			int h = left;
			for(size_t kernel_index = 0; kernel_index < kernel_width; ++kernel_index, ++h) {
				temp = i - h;
				if(temp < temp_vector[0].size() && temp >= 0 ) {
					product = (temp_vector[temp][j] * kernel[kernel_index]);
					convolve_sum += product;
				}
			}
			output[i][j] = convolve_sum;
			convolve_sum = 0, temp = 0;
		}

	}

	return output;
}


void smoothRow(Image* input, Image* output, vector<unsigned int> mask) {

	const unsigned int width = mask.size();
	size_t temp;
	size_t kernel_sum = 0;
	size_t convolve_sum = 0;
	const int left = -1 * (width / 2);

	for(unsigned int i = 0; i < width; ++i) {
		kernel_sum += mask[i];
	}

	for(size_t i = 0; i < input->num_rows(); ++i) {

		for(size_t j = 0; j < input->num_columns(); ++j) {
			// glide mask over
			int k = left;
			for(size_t kernel = 0; kernel < width; ++kernel, ++k) {
				temp = j - k;
				if(temp < input->num_columns() && temp >= 0 ) {
					convolve_sum += (mask[kernel] * input->GetPixel(i, temp));
				}
				else {
					convolve_sum += (mask[kernel] * input->GetPixel(i, j));
				}

			}
			output->SetPixel(i, j, convolve_sum/kernel_sum);
			convolve_sum = 0, temp = 0;
		}

	}

}

void smoothCol(Image* input, Image* output, vector<unsigned int> mask) {

	const unsigned int width = mask.size();
	size_t temp;
	size_t kernel_sum = 0;
	size_t convolve_sum = 0;
	const int top = -1 * (width / 2);

	for(unsigned int i = 0; i < width; ++i) {
		kernel_sum += mask[i];
	}

	for(size_t j = 0; j < input->num_columns(); ++j) {

		for(size_t i = 0; i < input->num_rows(); ++i) {
			// glide mask over
			int k = top;
			for(size_t kernel = 0; kernel < width; ++kernel, ++k) {
				temp = i - k;
				if(temp < input->num_rows() && temp >= 0 ) {
					convolve_sum += (mask[kernel] * input->GetPixel(temp, j));
				}
				else {
					convolve_sum += (mask[kernel] * input->GetPixel(i, j));
				}
			}
			output->SetPixel(i, j, convolve_sum/kernel_sum);
			convolve_sum = 0, temp = 0;
		}

	}


}

vector<unsigned int> getGaussianMask(Image* input, double sigma) {

	const unsigned int width = round(5*sigma);
	const int left = -1 * (width/2);
	double smallest_entry = gaussian(0, sigma, left);
	vector<unsigned int> mask;

	unsigned int value;
	int mask_x = left;

	// creating integer 1D gaussian kernel
	for(unsigned int i = 0; i < width; ++i) {
		value = round(gaussian(0, sigma, mask_x++) / smallest_entry);
		mask.push_back(value);
	}
	return mask;
}


Image findEdge(Image* input, double sigma, double tl, double th, double range) {

	vector<unsigned int> mask = getGaussianMask(input, sigma);

	vector<vector<int>> dx = {{-1, 1},{-1, 1}};
	vector<vector<int>> dy = {{1, 1},{-1, -1}};
	vector<vector<int>> gaussian_dx;
	vector<vector<int>> gaussian_dy;

	gaussian_dx = deriveGaussian(dx, mask);
	gaussian_dy = deriveGaussian(dy, mask);
/*
	for(const auto& i: gaussian_dx) {
		for(const auto& j: i) {
			cout << j << ' ';
		}
		cout << endl;
	}
*/
	long mask_sum = 0;
	long mask_width = mask.size();

	for(unsigned int i = 0; i < mask_width; ++i) {
		mask_sum += mask[i];
	}

	const double normalize_constant = 1.0/ ((mask_sum*mask_sum)*2.0);
	vector<vector<double>> dx_array = getFilteredArray(input, gaussian_dx, normalize_constant);
	vector<vector<double>> dy_array  = getFilteredArray(input, gaussian_dy, normalize_constant);


	//const double max_strength = abs(2*gaussian_dy[0][0]*255*normalize_constant);

	vector<vector<double>> strength(dx_array.size(), vector<double>(dx_array[0].size()));
	vector<vector<double>> orientation(dx_array.size(), vector<double>(dx_array[0].size()));

	Image edge_image;
	if(input->num_rows() == 0 || input->num_columns() == 0) {
		cout << "ERROR: empty image!" << endl;
		return edge_image;

	}

	// edge image is one row and one column smaller than actual image!!!
  edge_image.allocateFlatGrayScale(dx_array.size(), dx_array[0].size(), 0);
  edge_image.SetNumberGrayLevels(input->num_gray_levels());

	double dx_squared, dy_squared;


	double max_strength = 0;

	for(size_t i = 0; i < edge_image.num_rows(); ++i) {
		for(size_t j = 0; j < edge_image.num_columns(); ++j) {
			dx_squared = (dx_array[i][j] * dx_array[i][j]);
			dy_squared = (dy_array[i][j] * dy_array[i][j]);
			//printf("(dx, dy): %f, %f\n", dx_array[i][j], dy_array[i][j]);
			strength[i][j] = sqrt(dx_squared + dy_squared);
			if(strength[i][j] > max_strength)
				max_strength = strength[i][j];

			//printf("strength[%zu][%zu]: %f\n", i, j, strength[i][j]);
			orientation[i][j] = atan2(dy_array[i][j], dx_array[i][j]);
		}

	}

	nonmax_suppress(strength, orientation);


	// TODO: implement hysteresis threshold
/*
	double mapped_tl = (tl/range) * max_strength;
	double mapped_th = (th/range) * max_strength;


	set<list<pair<size_t, size_t>>> edge_set = hysteresis_threshold(strength,
																																	orientation,
																																	mapped_tl,
																																	mapped_th);

	for(const auto& i: edge_set) {
		for(const auto& j: i) {
			edge_image.SetPixel(j.first,
													j.second,
													round( (strength[j.first][j.second]/max_strength)*255) );

		}

	}
*/

	for(size_t i = 0; i < edge_image.num_rows(); ++i) {
		for(size_t j = 0; j < edge_image.num_columns(); ++j) {
			edge_image.SetPixel(i,
													j,
													round( (strength[i][j]/max_strength)*255) );
		}

	}

	return edge_image;

}

vector<vector<double>> getFilteredArray(Image* input, vector<vector<int>> kernel,
																								 			const double normal_constant) {

	// edge image is one row and one columns smaller than the actual image!!!
	vector<vector<double>> output(input->num_rows()-1, vector<double>(input->num_columns()-1));

	if(kernel.empty()) {
		cout << "Error: kernel is empty!" << endl;
		return output;
	}

	long convolve_sum = 0;

	for(size_t i = 0; i < input->num_rows()-1; ++i) {
		for(size_t j = 0; j < input->num_columns()-1; ++j) {

			for(size_t k_row = 0; k_row < kernel.size(); ++k_row) {
				for(size_t k_col = 0; k_col < kernel[0].size(); ++k_col) {

					convolve_sum += kernel[k_row][k_col] * input->GetPixel(i + k_row, j + k_col);
				}
			}
			output[i][j] = (convolve_sum*normal_constant);
			convolve_sum = 0;

		}

	}

	return output;

}

double fmod_math(double a, double b) {

	if( a < 0) {
		size_t temp = (a * -1) / b; // truncation
		return fmod( ((temp * b) + b + a), b);
	}
	return fmod(a, b);

}

unsigned short getClosestAngle(double angle) {
	// should return approximate angle that is normal to the given angle
	angle = fmod_math(angle, M_PI);

	if(angle < 0) {
		angle = angle + M_PI;
	}
	vector<double> diff_array = {abs(0 - angle),
												 			 abs((M_PI/4) - angle),
												 			 abs((M_PI/2) - angle),
												 			 abs((M_PI - (M_PI/4)) - angle)};

	double result = diff_array[0];
	unsigned short index = 0;

	for(size_t i = 1; i < diff_array.size(); ++i) {
		if(diff_array[i] < result) {
			result = diff_array[i];
			index = i;
		}
	}
	switch(index) {
		case 0: {
			return 0;}
		case 1: {
			return 45;}
		case 2: {
			return 90;}
		case 3: {
			return 135;}
	}

}


void nonmax_suppress(vector<vector<double>>& strength,
										 const vector<vector<double>>& orientation) {

	vector<vector<double>> temp(strength);
	vector<vector<unsigned short>> angle_array(strength.size(),
																						 vector<unsigned short>(strength[0].size()));

	unsigned short angle;

	for(size_t i = 0; i < strength.size(); ++i) {
		for(size_t j = 0; j < strength[0].size(); ++j) {
			angle = getClosestAngle(orientation[i][j]);
			switch(angle) {
				case 0: {

					// if bordering pixel is inside image
					if(j+1 < strength[0].size() && j != 0) {
						if( (strength[i][j] < strength[i][j+1]) ||
									(strength[i][j] < strength[i][j-1]) ) {
							temp[i][j] = 0;
						}

					}
					else{ // a bordering pixel is outside of image
						if( (j+1 >= strength[0].size()) &&
								(j != 0) &&
								(strength[i][j] < strength[i][j-1]) ) {
							temp[i][j] = 0;
						}
					 else if( (j == 0) &&
										(j+1 < strength[0].size()) &&
										(strength[i][j] < strength[i][j+1]) ) {
							temp[i][j] = 0;
						}
					}

					break;
				}
				case 45: {
					// if current pixel neighbors are inside the image
					if( (i != 0 && j+1 < strength[0].size()) &&
							(i+1 < strength.size() && j != 0)) {

						// if pixel strength is less than any of its diagonal neighbors
						if( (strength[i][j] < strength[i-1][j+1]) ||
								(strength[i][j] < strength[i+1][j-1]) ) {
							temp[i][j] = 0;
						}

					}
					else{ // not at bottom right corner
						// a bordering pixel is out of the image
						// top-right pixel is out of image
						if( (i == 0 || j+1 >= strength[0].size()) &&
								(i+1 < strength.size()) &&
								(j != 0) &&
								(strength[i][j] < strength[i+1][j-1]) ) {
							temp[i][j] = 0;
						}
						else if( (i+1 >= strength.size() || j == 0) &&
										 (i != 0) &&
										 (j+1 < strength[0].size()) &&
										 (strength[i][j] < strength[i-1][j+1]) ) {
							// bottom-left pixel is out of image
							temp[i][j] = 0;
						}

					}
					break;
				}
				case 90: {
					// if bordering pixel is out of image
					if(i+1 < strength.size() && i != 0) {
						if( (strength[i][j] < strength[i+1][j]) ||
									(strength[i][j] < strength[i-1][j]) ) {
							temp[i][j] = 0;
						}
					}
					else { // a bordering pixel is at the edge of image
						if( (i+1 >= strength.size()) &&
								 (i != 0) &&
								(strength[i][j] < strength[i-1][j]) ) {
							temp[i][j] = 0;
						}
					 	else if( (i == 0) &&
										(i+1 < strength.size()) &&
										(strength[i][j] < strength[i+1][j]) ) {
							temp[i][j] = 0;
						}
					}

					break;
				}
				case 135: {
					// if no bordering pixels are outside of image
					if( (i != 0 && j != 0) &&
							(i+1 < strength.size() && j+1 < strength[0].size())) {
						// if pixel strength is less than any of its diagonal neighbors

						if( (strength[i][j] < strength[i-1][j-1]) ||
								(strength[i][j] < strength[i+1][j+1]) ) {

							temp[i][j] = 0;
						}
					}
					else{
						// a bordering pixel is out of the image
						// top-left pixel is out of image
						if( (i == 0 || j == 0) &&
								(i+1 < strength.size()) &&
								(j+1 < strength[0].size()) &&
								(strength[i][j] < strength[i+1][j+1]) ) {

							temp[i][j] = 0;
						}
						else if( (i+1 >= strength.size() || j+1 >= strength[0].size()) &&
										 (i != 0) &&
										 (j != 0) &&
										 (strength[i][j] < strength[i-1][j-1]) ) {
							// bottom right pixel is out of image
							temp[i][j] = 0;
						}

					}
					break;
				} //case 135

			}
		}

	}


	strength = temp;

}
// TODO: implement hysteresis thresholding
/*
set<list<pair<size_t, size_t>>> hysteresis_threshold(vector<vector<double>>& strength,
													const vector<vector<double>>& orientation,
													double tl, double th) {


	set<list<pair<size_t, size_t>>> edge_set;
	list<pair<size_t, size_t>> edge_list;
	pair<size_t, size_t> coordinate;

	if(tl > th) {
		cout << "Invalid threshold values! Condition must be: (tl < th)" << endl;
		cout << "Returning empty set!" << endl;
		return edge_set;
	}

	vector<vector<bool>> visited(orientation.size(),
															 vector<bool>(orientation[0].size()));

	vector<vector<double>> edge_orientation(orientation.size(),
																					vector<double>(orientation[0].size(), false));

	size_t x, y;
	bool top_left, top, top_right, left, right, bottom_left, bottom, bottom_right;
	pair<size_t, size_t> path_1;
	pair<size_t, size_t> path_2;
	short window[] = {-1, 0, 1};

	double current_strength = 0.0;

	for(size_t i = 0; i < strength.size(); ++i) {
		for(size_t j = 0; j < strength.size(); ++j) {

			if(strength[i][j] > th && !visited[i][j]) {
				visited[i][j] = true;

				// look for neighbors
				for(unsigned short h = 0; h < 3; ++h) {
					for(unsigned short k = 0; k < 3; ++k) {
						if(i == 0 && j == 0) {


						}
						else if(i == (strength.size()-1) && j == (strength[0].size()-1)) {


						}
						else if(i == 0) {


						}
						else if(j == 0) {

						}
						else if(i == (strength.size()-1) ) {


						}
						else if(j == (strength[0].size()-1)) {

						}

					}

				}


				while(){
					// top row (i-1)
					if(i != 0 && j !=0 && strength[i-1][j-1] > tl && strength[i-1][j-1] > current_strength &&
																																								!visited[i-1][j-1]) {
						current_path.first = i-1;
						current_path.second = j-1;
						current_strength = strength[i-1][j-1];
						visited[i-1][j-1] = true;
					}
					if(i != 0 && strength[i-1][j] > tl && strength[i-1][j] > current_strength) {
						current_path.first = i-1;
						current_path.second = j;
						current_strength = strength[i-1][j-1];
						visited[i-1][j] = true;
					}
					if(i != 0 && (j+1) < strength[0].size() && strength[i-1][j+1] > tl &&
																										 strength[i-1][j+1] > current_strength &&
																										 !visited[i-1][j+1]) {
						current_path.first = i-1;
						current_path.second = j+1;
						current_strength = strength[i-1][j+1];
						visited[i-1][j+1] = true;

					}
					// middle row (same row as current)
					if(j != 0 && strength[i][j-1] > tl && strength[i][j-1] > current_strength &&
																								!visited[i][j-1]) {
						current_path.first = i;
						current_path.second = j-1;
						current_strength = strength[i][j-1];
						visited[i][j-1] = true;
					}
					if((j+1) < strength[0].size() && strength[i][j+1] > tl &&
																					 strength[i][j+1] > current_strength &&
																					 !visited[i][j+1]	) {
						current_path.first = i;
						current_path.second = j+1;
						current_strength = strength[i][j+1];
						visited[i][j+1] = true;
					}
					// bottom row (i = i+1)
					if((i+1) < strength.size() && j != 0 && strength[i+1][j-1] > tl &&
																									strength[i+1][j-1] > current_strength &&
																									!visited[i+1][j-1]) {
						current_path.first = i+1;
						current_path.second = j-1;
						current_strength = strength[i+1][j-1];
						visited[i+1][j-1] = true;

					}
					if((i+1) < strength.size() && strength[i+1][j] > tl && strength[i+1][j] > current_strength &&
																																 !visited[i+1][j]) {
						current_path.first = i+1;
						current_path.second = j;
						current_strength = strength[i+1][j];
						visited[i+1][j] = true;


					}
					if((i+1) < strength.size() && (j+1) < strength[0].size() &&
																								strength[i+1][j+1] > tl &&
																								strength[i+1][j+1] > current_strength &&
																								!visited[i+1][j+1]) {
						current_path.first = i+1;
						current_path.second = j+1;
						current_strength = strength[i+1][j+1];
						visited[i+1][j+1] = true;
					}


					list.push_back(current_path);

				} // end while

				edge_set.insert(edge_list);
				edge_list = list<pair<size_t, size_t>>();
			}

		}

	}


}
*/
void getHoughSpace(Image* edge_image, Image* hough_image,
																					 const size_t theta_resolution,
																					 const size_t rho_resolution,
																					 const string hough_array_file) {

	// xcost + ysint = p
	// f(t) = p = xcost + ysint
	// t is in [0, pi]
	const double image_diagonal = sqrt( (edge_image->num_rows() * edge_image->num_rows()) +
																		 	(edge_image->num_columns() * edge_image->num_columns()) );

	const double theta_increment = M_PI / theta_resolution;
	const double rho_increment = image_diagonal / (rho_resolution / 2);
	//pair<double, double> theta_rho_increment_pair(theta_increment, rho_increment);
	// IMPORTANT: rho can be negative, so rho_resolution is divided
	// by 2, half above theta axis and half below it
	const size_t theta_axis = round(rho_resolution / 2);

	size_t discrete_theta = 0;
	long discrete_rho;
	double theta = 0, rho;
	unsigned short image_vote;

	vector<vector<size_t>> hough_array(rho_resolution, vector<size_t>(theta_resolution, 0));

	size_t largest_vote = 0;

	for(size_t i = 0; i < edge_image->num_rows(); ++i) {
		for(size_t j = 0; j < edge_image->num_columns(); ++j) {
			// check if pixel is an edge
			if(edge_image->GetPixel(i, j) == 255) {
				// x = i, y = j
				while(discrete_theta < theta_resolution) {
					rho = i*cos(theta) + j*sin(theta);
					discrete_rho = theta_axis - round(rho / rho_increment);

					if(discrete_rho < rho_resolution && discrete_rho > 0) {
						hough_array[discrete_rho][discrete_theta]++;

						if(hough_array[discrete_rho][discrete_theta] > largest_vote)
							largest_vote = hough_array[discrete_rho][discrete_theta];
					}
					theta += theta_increment;
					discrete_theta++;
				}
				discrete_theta = 0;
				theta = 0;

			} // end if
		}
	}

	unsigned short relative_gray_value;

	// initialize hough image
	hough_image->allocateFlatGrayScale(rho_resolution, theta_resolution, 0);
	hough_image->SetNumberGrayLevels(edge_image->num_gray_levels());
	ofstream data;
	data.open(hough_array_file);

	if(data.is_open()) {
		data.precision(12);
		data << theta_axis << ' ' << theta_resolution << ' ' << rho_resolution << ' '
				 << theta_increment << ' ' << rho_increment << '\n';

		for(size_t i = 0; i < hough_array.size(); ++i) {
			for(size_t j = 0; j < hough_array[0].size(); ++j) {
				// write to file
				data << hough_array[i][j] << ' ';
				// write to image
				relative_gray_value = round(((double)hough_array[i][j] / (double)largest_vote)*255) ;
				hough_image->SetPixel(i, j, relative_gray_value);

			}
			data << '\n';
		}
	}
	else {
		cout << "ERROR: can't open file " << hough_array_file << '!' << endl;
		cout << hough_array_file << " is unchanged!" << endl;
		//throw "ERROR: Unable to open " << hough_array_file  << '!' << endl;
		return;
	}

	data.close();

}

Image getHoughImage(const string& hough_array_file, const unsigned short graylevel) {

	ifstream data;
	data.open(hough_array_file);

	Image hough_image;
	hough_image.SetNumberGrayLevels(graylevel);

	if(data.is_open()) {
		// collect initial data
		// <theta resolution, rho resolution, theta increment, rho increment>
		size_t theta_axis, theta_resolution, rho_resolution;
		double theta_increment, rho_increment;
		data.precision(12);
		data >> theta_axis;
		data.ignore(256, ' ');
		data >> theta_resolution;
		data.ignore(256, ' ');
		data >> rho_resolution;
		data.ignore(256, ' ');
		data >> theta_increment;
		data.ignore(256, ' ');
		data >> rho_increment;
		data.ignore(256, ' ');

		//printf("%zu %zu %.12f %.12f\n", theta_resolution, rho_resolution, theta_increment, rho_increment);

		vector<vector<size_t>> hough_array(rho_resolution, vector<size_t>(theta_resolution));
		size_t largest_vote = 0;

		for(size_t i = 0; i < rho_resolution; ++i) {
			for(size_t j = 0; j < theta_resolution; ++j) {
				data >> hough_array[i][j];
				if(hough_array[i][j] > largest_vote) {
					largest_vote = hough_array[i][j];
				}
				data.ignore(256, ' ');
			}
		}

		hough_image.allocateFlatGrayScale(rho_resolution, theta_resolution, 0);
		unsigned short relative_gray_value;

		for(size_t i = 0; i < hough_array.size(); ++i) {
			for(size_t j = 0; j < hough_array[0].size(); ++j) {
				relative_gray_value = round(((double)hough_array[i][j] / (double)largest_vote)*255);
				hough_image.SetPixel(i, j, relative_gray_value);

			}
		}
	}
	else {
		cout << "ERROR: can't open file " << hough_array_file << '!' << endl;
		cout << "The output Hough image is empty!" << endl;
	}

	return hough_image;
}

void getHoughPatchSpace(Image* edge_image, Image* hough_image,
																					 const size_t theta_resolution,
																					 const size_t rho_resolution,
																					 const string hough_array_file) {

	const double image_diagonal = sqrt( (edge_image->num_rows() * edge_image->num_rows()) +
																		 	(edge_image->num_columns() * edge_image->num_columns()) );

	const double theta_increment = M_PI / theta_resolution;
	const double rho_increment = image_diagonal / (rho_resolution / 2);

	const size_t theta_axis = round(rho_resolution / 2);

	size_t discrete_theta = 0;
	long discrete_rho;
	double theta = 0, rho;
	unsigned short image_vote;

	vector<vector<size_t>> hough_array(rho_resolution, vector<size_t>(theta_resolution, 0));

	size_t largest_vote = 0;

	for(size_t i = 0; i < edge_image->num_rows(); ++i) {
		for(size_t j = 0; j < edge_image->num_columns(); ++j) {
			// check if pixel is an edge
			if(edge_image->GetPixel(i, j) == 255) {
				// x = i, y = j
				while(discrete_theta < theta_resolution) {
					rho = i*cos(theta) + j*sin(theta);

					// rounding error carries over to line drawing
					discrete_rho = theta_axis - round(rho / rho_increment);

					if(discrete_rho < rho_resolution && discrete_rho > 0) {
						hough_array[discrete_rho][discrete_theta]++;

						// get largest vote
						if(hough_array[discrete_rho][discrete_theta] > largest_vote)
							largest_vote = hough_array[discrete_rho][discrete_theta];

						// upper row
						if(discrete_rho != 0 && discrete_theta != 0) {
							hough_array[discrete_rho-1][discrete_theta-1]++;

							if(hough_array[discrete_rho-1][discrete_theta-1] > largest_vote)
								largest_vote = hough_array[discrete_rho-1][discrete_theta-1];
						}
						if(discrete_rho != 0) {
							hough_array[discrete_rho-1][discrete_theta]++;

							if(hough_array[discrete_rho-1][discrete_theta] > largest_vote)
								largest_vote = hough_array[discrete_rho-1][discrete_theta];
						}
						if(discrete_rho != 0 && discrete_theta < (theta_resolution-1)) {
							hough_array[discrete_rho-1][discrete_theta+1]++;

							if(hough_array[discrete_rho-1][discrete_theta+1] > largest_vote)
								largest_vote = hough_array[discrete_rho-1][discrete_theta+1];
						}

						// middle row
						if(discrete_theta != 0) {
							hough_array[discrete_rho][discrete_theta-1]++;

							if(hough_array[discrete_rho][discrete_theta-1] > largest_vote)
								largest_vote = hough_array[discrete_rho][discrete_theta-1];
						}
						if(discrete_theta < (theta_resolution+1)) {
							hough_array[discrete_rho][discrete_theta+1]++;

							if(hough_array[discrete_rho][discrete_theta+1] > largest_vote)
								largest_vote = hough_array[discrete_rho][discrete_theta+1];
						}

						// lower row
						if(discrete_rho < (rho_resolution - 1) && discrete_theta != 0) {
							hough_array[discrete_rho+1][discrete_theta-1]++;

							if(hough_array[discrete_rho+1][discrete_theta-1] > largest_vote)
								largest_vote = hough_array[discrete_rho+1][discrete_theta-1];

						}
						if(discrete_rho < (rho_resolution - 1)) {
							hough_array[discrete_rho+1][discrete_theta]++;

							if(hough_array[discrete_rho+1][discrete_theta] > largest_vote)
								largest_vote = hough_array[discrete_rho+1][discrete_theta];
						}
						if(discrete_rho < (rho_resolution - 1) && discrete_theta < (theta_resolution - 1)) {
							hough_array[discrete_rho+1][discrete_theta+1]++;

							if(hough_array[discrete_rho+1][discrete_theta+1] > largest_vote)
								largest_vote = hough_array[discrete_rho+1][discrete_theta+1];
						}


					}
					theta += theta_increment;
					discrete_theta++;
				}
				discrete_theta = 0;
				theta = 0;

			} // end if

		}

	}

	unsigned short relative_gray_value;

	// initialize hough image
	hough_image->allocateFlatGrayScale(rho_resolution, theta_resolution, 0);
	hough_image->SetNumberGrayLevels(edge_image->num_gray_levels());
	ofstream data;
	data.open(hough_array_file);

	if(data.is_open()) {
		data.precision(12);
		data << theta_axis << ' ' << theta_resolution << ' ' << rho_resolution << ' '
				 << theta_increment << ' ' << rho_increment << '\n';

		for(size_t i = 0; i < hough_array.size(); ++i) {
			for(size_t j = 0; j < hough_array[0].size(); ++j) {
				// write to file
				data << hough_array[i][j] << ' ';
				// write to image
				relative_gray_value = round(((double)hough_array[i][j] / (double)largest_vote)*255) ;
				hough_image->SetPixel(i, j, relative_gray_value);

			}
			data << '\n';
		}
	}
	else {
		cout << "ERROR: can't open file " << hough_array_file << '!' << endl;
		cout << hough_array_file << " is unchanged!" << endl;
		//throw "ERROR: Unable to open " << hough_array_file  << '!' << endl;
		return;
	}

	data.close();

}

unordered_map<unsigned short, double> getNonhomogenousCenterX(Image* labeled_hough,
																															Image* original_hough) {

	unordered_map<unsigned short, unsigned long> label_to_area = getNonhomogenousArea(labeled_hough,
																																										original_hough);
	unordered_map<unsigned short, double> label_to_x_moment;
	unordered_map<unsigned short, double>::const_iterator it;

	int pixel;

	for(size_t i = 0; i < labeled_hough->num_rows(); ++i) {

		for(size_t j = 0; j < labeled_hough->num_columns(); ++j) {

			pixel = labeled_hough->GetPixel(i, j);
			it = label_to_x_moment.find(pixel);

			if(pixel != 0) {

				if(it == label_to_x_moment.end() ) {
					label_to_x_moment[pixel] = i*original_hough->GetPixel(i, j);
				}
				else {
					label_to_x_moment[pixel] += i*original_hough->GetPixel(i, j);

				}

			}


		}
	}

	unsigned long area;
	for(auto& pair: label_to_x_moment) {
		area = label_to_area[pair.first];
		pair.second /= area;

	}

	return label_to_x_moment;

}

unordered_map<unsigned short, double> getNonhomogenousCenterY(Image* labeled_hough,
																															Image* original_hough) {

	unordered_map<unsigned short, unsigned long> label_to_area = getNonhomogenousArea(labeled_hough,
																																										original_hough);
	unordered_map<unsigned short, double> label_to_y_moment;
	unordered_map<unsigned short, double>::const_iterator it;

	int pixel;

	for(size_t i = 0; i < labeled_hough->num_rows(); ++i) {

		for(size_t j = 0; j < labeled_hough->num_columns(); ++j) {

			pixel = labeled_hough->GetPixel(i, j);
			it = label_to_y_moment.find(pixel);

			if(pixel != 0) {

				if(it == label_to_y_moment.end() ) {
					label_to_y_moment[pixel] = j*original_hough->GetPixel(i, j);
				}
				else {
					label_to_y_moment[pixel] += j*original_hough->GetPixel(i, j);

				}

			}


		}
	}

	unsigned long area;
	for(auto& pair: label_to_y_moment) {
		area = label_to_area[pair.first];
		pair.second /= area;

	}

	return label_to_y_moment;

}

list<pair<double, double>> getLineParameters(const unsigned short threshold,
																						 const string hough_array_file,
																						 const double theta_error,
																						 const double rho_error,
																						 const unsigned short graylevel,
																						 const unsigned short base,
																						 const unsigned short increment) {
	ifstream data;
	data.open(hough_array_file);
	list<pair<double, double>> theta_rho_pairs;

	if(data.is_open()) {

		Image original_hough, labeled_hough_image;

		size_t theta_axis, theta_resolution, rho_resolution;
		double theta_increment, rho_increment;
		data.precision(12);
		data >> theta_axis;
		data.ignore(256, ' ');
		data >> theta_resolution;
		data.ignore(256, ' ');
		data >> rho_resolution;
		data.ignore(256, ' ');
		data >> theta_increment;
		data.ignore(256, ' ');
		data >> rho_increment;
		data.ignore(256, ' ');

		original_hough = getHoughImage(hough_array_file, graylevel);
		labeled_hough_image = original_hough;
		toBinary(&labeled_hough_image, threshold);

		try {
			separateObjectsByShade(&labeled_hough_image, base, increment);
		}
		catch(const char* msg) {
			cerr << msg << endl;

		}

		// better algorithm:
		// use a labeled image alongside original hough image
		// to get central moment of non-homogenous body
		unordered_map<unsigned short, double> nonhomogenousCenterX = getNonhomogenousCenterX(
																																 	&labeled_hough_image,
																																	&original_hough);

		unordered_map<unsigned short, double> nonhomogenousCenterY = getNonhomogenousCenterY(
																																 	&labeled_hough_image,
																																	&original_hough);
		//labeled_hough_image.update();
		//unordered_map<unsigned short, double> centerX = labeled_hough_image.getCenterX();
		//unordered_map<unsigned short, double> centerY = labeled_hough_image.getCenterY();

		// rho is X, theta is Y
		double x, y;
		double approximate_rho;

		for(const auto& cur: nonhomogenousCenterX) {
			x = cur.second;
			y = nonhomogenousCenterY[cur.first];

			approximate_rho = (theta_axis - x) * rho_increment;

			theta_rho_pairs.push_back(pair<double, double>(y*theta_increment, approximate_rho));
			//printf("(theta, rho): %f, %f\n", y*theta_increment, approximate_rho);
			//printf("(x, y): %f, %f\n", y, x);
			//printf("(x, y): %f, %f\n\n", y, theta_axis - x);
		}

		theta_rho_pairs = filterRedundantLines(&original_hough, theta_rho_pairs, theta_axis, theta_increment, rho_increment, theta_error, rho_error);
		//(*labeled_hough) = labeled_hough_image;
		data.close();
	}

	else {
		cout << "ERROR: can't open file " << hough_array_file << '!' << endl;
		cout << "The output Hough image is empty!" << endl;

	}


	return theta_rho_pairs;

}
// TODO: Not necessary and miniizing might reduce information on a complex image
list<pair<double, double>> filterRedundantLines(Image* hough, list<pair<double, double>> theta_rho_pairs,
																								size_t theta_axis, double theta_increment, double rho_increment,
																								unsigned int theta_inc_multiple_error,
																								unsigned int rho_inc_multiple_error) {


	list<pair<double, double>> filtered_parameters = theta_rho_pairs;

	list<pair<double, double>>::iterator outer_it = filtered_parameters.begin();
	list<pair<double, double>>::iterator inner_it = filtered_parameters.begin();
	list<pair<double, double>>::iterator temp;

	size_t x, y, x2, y2, param_size = filtered_parameters.size();
	double theta_error, rho_error;

	double theta, rho, theta2, rho2, best_theta, best_rho;

	while(outer_it != filtered_parameters.end()) {
		theta = outer_it->first;
		rho = outer_it->second;
		if(theta > M_PI-0.05 && rho < 0){

			outer_it = filtered_parameters.erase(outer_it);
			continue;

		}
		else if(theta > M_PI-0.05 && rho >= 0) {
			outer_it->first = 0.0;
		}
		outer_it++;

	}

	outer_it = filtered_parameters.begin();

	for(; outer_it != filtered_parameters.end();) {
		bool found_better_pair = false;

		theta = outer_it->first;
		rho = outer_it->second;

		inner_it = outer_it;
		inner_it++;

		best_theta = theta;
		best_rho = rho;

		// get (x, y), the center of hough section
		x = theta_axis - (rho/rho_increment);
		y = theta / theta_increment;

		for(; inner_it != filtered_parameters.end(); ++inner_it) {
			if(inner_it == outer_it) {
				continue;
			}
			theta2 = inner_it->first;
			rho2 = inner_it->second;

			theta_error = abs(best_theta - theta2);
			rho_error = abs(best_rho - rho2);
/*
			printf("theta: %f, rho: %f\n", theta*(180.0/M_PI), rho);
			printf("theta2: %f, rho2: %f\n", theta2*(180.0/M_PI), rho2);
			printf("theta_error: %f, rho_error: %f\n", theta_error, rho_error);
			printf("theta bound: %f, rho bound: %f\n", theta_inc_multiple_error*(theta_increment), rho_inc_multiple_error*(rho_increment));
*/
			// if within error bounds, find best line
			if(theta_error <= theta_inc_multiple_error*(theta_increment) &&
				 rho_error <= rho_inc_multiple_error*(rho_increment)) {
				// get another (x, y), the center of another hough section
				x2 = theta_axis - (rho2/rho_increment);
				y2 = theta2 / theta_increment;

				if(hough->GetPixel(x,y) == hough->GetPixel(x2,y2)) {

					found_better_pair = true;
					outer_it = filtered_parameters.erase(outer_it);

				}
				if(hough->GetPixel(x,y) < hough->GetPixel(x2,y2)) {

					//temp = outer_it;
					found_better_pair = true;
					outer_it = filtered_parameters.erase(outer_it);
					//param_size--;
				}
				else if(hough->GetPixel(x,y) > hough->GetPixel(x2,y2)) {
					//temp = outer_it;

					found_better_pair = true;
					inner_it = filtered_parameters.erase(inner_it);
					//param_size--;
				}

			}
		//cout << endl;
		}

		if(!found_better_pair)
			++outer_it;

	}

	return filtered_parameters;

}

Image getSectionedHough(const unsigned short threshold,
										  const string hough_array_file,
										  const unsigned short graylevel,
										  const unsigned short base,
										  const unsigned short increment) {
	ifstream data;
	data.open(hough_array_file);

	Image output_hough, labeled_hough;

	if(data.is_open()) {

		size_t theta_axis, theta_resolution, rho_resolution;
		double theta_increment, rho_increment;
		data.precision(12);
		data >> theta_axis;
		data.ignore(256, ' ');
		data >> theta_resolution;
		data.ignore(256, ' ');
		data >> rho_resolution;
		data.ignore(256, ' ');
		data >> theta_increment;
		data.ignore(256, ' ');
		data >> rho_increment;
		data.ignore(256, ' ');

		output_hough = getHoughImage(hough_array_file, graylevel);
		labeled_hough = output_hough;
		toBinary(&labeled_hough, threshold);

		try {
			separateObjectsByShade(&labeled_hough, base, increment);
		}
		catch(const char* msg) {
			cerr << msg << endl;

		}

		// rho is X, theta is Y
		double x, y;
		double approximate_rho;

		for(size_t i = 0; i < output_hough.num_rows(); ++i) {
			for(size_t j = 0; j < output_hough.num_columns(); ++j) {
				if(labeled_hough.GetPixel(i, j) == 0)
					output_hough.SetPixel(i, j, 0);
			}


		}
		data.close();
	}

	else {
		cout << "ERROR: can't open file " << hough_array_file << '!' << endl;
		cout << "The output Hough image is empty!" << endl;

	}


	return output_hough;

}

void drawDetectedLines(Image* input, list<pair<double, double>> theta_rho_pairs) {

	size_t height = input->num_rows();
	size_t width = input->num_columns();
	double theta, rho, rounded_rho;
	size_t x, y;
	double constant, coefficient;

	// this is basically a n * 4 array
	// each vector contains two points {x1, y1, x2, y2}
	// (x, y) = (i, j) where i is row number and j is column number
	list<list<size_t>> point_pairs_list;
	list<size_t> point_pair;

	for(const auto& i: theta_rho_pairs) {
		theta = i.first;
		rho = i.second;


		// vertical line when normal vector to the line is 0 radian
		if(theta < 0.01 && theta >= 0.0 || theta <= M_PI && theta > (M_PI - 0.01)) {
		//if(theta <  0.05 && theta >= 0.0) {
			// x = (rho / cos(theta)) - y*( sin(theta)/cos(theta) )
			// let y = 0, and since theta ~ 0
			// so cos(theta) ~ 1, therefore x ~ rho
			if(rho < 0 || rho >= height) // line is outside of image
				continue;

			//printf("1. theta: %f, (%f, %f), (%f, %zu)\n", theta, rho, 0.0, rho, (width-1));
			rounded_rho = round(rho);

			point_pair.push_back(rounded_rho);
			point_pair.push_back(0);

			point_pair.push_back(rounded_rho);
			point_pair.push_back(width-1);

		}
		// horizontal line when normal vector to the line is pi/2 radians
		else if(theta < ((M_PI/2) + 0.01) && theta > ((M_PI/2)-0.01)) {
		//if(theta <  0.05 && theta >= 0.0) {
			// x = (rho / cos(theta)) - y*( sin(theta)/cos(theta) )
			// let y = 0, and since theta ~ 0
			// so cos(theta) ~ 1, therefore x ~ rho


			rounded_rho = round(rho);

			point_pair.push_back(0);
			point_pair.push_back(rounded_rho);

			point_pair.push_back(height-1);
			point_pair.push_back(rounded_rho);

		}

		else {
			// y = (rho / sin(t)) - x*(cos(t)/sin(t))
/*
			printf("2. theta: %f, (%f, %f), (%zu, %f)\n", theta, 0.0, rho / sin(theta),
																				(height-1),
																				(rho / sin(theta)) - (height-1)*(cos(theta)/sin(theta)));
*/

			bool found_first_point = false;

			x = 0;
			y = round(rho / sin(theta));

			if(y >= width || y < 0) {
				x++;
				y = round((rho / sin(theta))) - x*(cos(theta)/sin(theta));

				while((y >= width || y < 0) && x < height) {
					x++;
					y = round((rho / sin(theta))) - x*(cos(theta)/sin(theta));
				}

			}
			//printf("(x, y): %zu, %zu\n", x, y);
			if(x < height) {
				found_first_point = true;
				point_pair.push_back(x);
				point_pair.push_back(y);
			}

			x = height-1;
			y = round((rho / sin(theta))) - x*(cos(theta)/sin(theta));

			if(y >= width || y < 0) {
				x--;
				y = round(rho / sin(theta)) - x*(cos(theta)/sin(theta));
				while((y >= width || y < 0) && x >= 0) {
					x--;
					y = round((rho / sin(theta))) - x*(cos(theta)/sin(theta));
				}

			}
			//printf("(x, y): %zu, %zu\n\n", x, y);
			if(found_first_point) {
				point_pair.push_back(x);
				point_pair.push_back(y);
			}


		}

		point_pairs_list.push_back(point_pair);
		point_pair = list<size_t>();
	}



	list<size_t>::const_iterator it;
	size_t x1, y1, x2, y2;

	for(const auto& i: point_pairs_list) {
		it = i.cbegin();
		x1 = *it++;
		y1 = *it++;
		x2 = *it++;
		y2 = *it;

		//printf("(x1, y1), (x2, y2): (%zu, %zu), (%zu, %zu)\n\n", x1, y1, x2, y2);
		DrawLine(x1, y1, x2, y2, 255, input);

	}



}

// TODO: get draw lines that match the object edges
/*
void drawEdgeLines(Image* input, const Image* binary_image, list<pair<double, double>> theta_rho_pairs){

	size_t height = input->num_rows();
	size_t width = input->num_columns();
	double theta, rho, rounded_rho;
	size_t x, y;
	double constant, coefficient;

	// this is basically a n * 4 array
	// each vector contains two points {x1, y1, x2, y2}
	// (x, y) = (i, j) where i is row number and j is column number
	list<list<size_t>> point_pairs_list;
	list<size_t> point_pair;
	bool starting_point_found = false;
	size_t first_x, first_y;

	for(const auto& i: theta_rho_pairs) {
		theta = i.first;
		rho = i.second;

		// vertical line when normal vector to the line is 0 radian
		if(theta <  0.03 && theta >= 0.0 || theta <=  M_PI && theta > (M_PI - 0.03)) {
			// x = (rho / cos(theta)) - y*( sin(theta)/cos(theta) )
			// let y = 0, and since theta ~ 0
			// so cos(theta) ~ 1, therefore x ~ rho
			if(rho < 0 || rho >= height) // line is outside of image
				continue;

			//printf("1. theta: %f, (%f, %f), (%f, %zu)\n", theta, rho, 0.0, rho, (width-1));
			rounded_rho = round(rho);

			first_x = rounded_rho;
			first_y = 0;


			// check to left and right of line orthogonally
			while()


			point_pair.push_back(rounded_rho);
			point_pair.push_back(0);

			point_pair.push_back(rounded_rho);
			point_pair.push_back(width-1);

		}

		while(!starting_point_found) {


		}

		point_pairs_list.push_back(point_pair);
		point_pair = list<size_t>();
	}



	list<size_t>::const_iterator it;
	size_t x1, y1, x2, y2;

	for(const auto& i: point_pairs_list) {
		it = i.cbegin();
		x1 = *it++;
		y1 = *it++;
		x2 = *it++;
		y2 = *it;

		//printf("(x1, y1), (x2, y2): (%zu, %zu), (%zu, %zu)\n\n", x1, y1, x2, y2);
		DrawLine(x1, y1, x2, y2, 255, input);

	}


}

*/

list<pair<double, double>> getIntersections(Image* image, list<pair<double, double>> theta_rho_pairs) {

	double theta, rho, theta2, rho2, numerator, denominator, x, y;
	list<pair<double, double>> intersections;
	list<pair<double, double>>::iterator out_it = theta_rho_pairs.begin();
	list<pair<double, double>>::iterator in_it;

	size_t parameter_list_size = theta_rho_pairs.size();
	// x = (rho / cos(t)) - y*(sin(t)/cos(t))
	// y = (rho / sin(t)) - x*(cos(t)/sin(t))

	size_t total_rows = image->num_rows();
	size_t total_columns = image->num_columns();


	cout << "theta rho pairs size: " << theta_rho_pairs.size() << endl;

	// iterate through theta_rho_pairs
	for(; out_it != theta_rho_pairs.end(); out_it = theta_rho_pairs.erase(out_it) ) {
		theta = out_it->first;
		rho = out_it->second;

		//printf("out_it(theta, rho): (%f, %f)\n", (out_it->first)*(180.0/M_PI), out_it->second);

		for(in_it = theta_rho_pairs.begin(); in_it != theta_rho_pairs.end(); ++in_it) { // less 1 since we skip current value
			if(in_it == out_it) {
				continue;
			}

			//printf("in_it(theta, rho): (%f, %f)\n", (in_it->first)*(180.0/M_PI), in_it->second);
			theta2 = in_it->first;
			rho2 = in_it->second;

			// when theta ~ 0 or theta ~ M_PI, sin(theta) ~ 0
			// so we want to use cost in denominator
			//if(theta <  0.05 && theta >= 0.0 || theta <=  M_PI && theta > M_PI-0.05) {
			if((theta <  0.05 && theta >= 0.0 || theta <=  M_PI && theta > M_PI-0.05) &&
					(theta2 >= 0.5 || theta2 <= M_PI-0.05)) {

				//printf("theta ~ 0, 180\n");
				numerator = (rho/cos(theta)) - (rho2/cos(theta2));
				denominator = (sin(theta)/cos(theta)) - (sin(theta2)/cos(theta2));
				y = numerator/denominator;
				x = (rho / cos(theta)) - y*(sin(theta)/cos(theta));

			}
			else if((theta2 <  0.05 && theta2 >= 0.0 || theta2 <=  M_PI && theta2 > M_PI-0.05) &&
					(theta >= 0.5 || theta <= M_PI-0.05)) {

				//printf("theta ~ 0, 180\n");
				numerator = (rho/cos(theta)) - (rho2/cos(theta2));
				denominator = (sin(theta)/cos(theta)) - (sin(theta2)/cos(theta2));
				y = numerator/denominator;
				x = (rho / cos(theta)) - y*(sin(theta)/cos(theta));

			}
			else {
				//printf("theta !~ 0, 180\n");
				numerator = (rho/sin(theta)) - (rho2/sin(theta2));
				denominator = (cos(theta)/sin(theta)) - (cos(theta2)/sin(theta2));
				x = numerator/denominator;
				y = (rho / sin(theta)) - x*(cos(theta)/sin(theta));
			}

			//printf("x: %f , y: %f\n", x, y);
			if(rho < 0 || x >= total_rows || x < 0 || y >= total_columns || y < 0 || denominator == 0) {
				continue;
			}

			//printf("num: %f , denom: %f\n", numerator, denominator);
			//printf("chosen(x, y): (%f, %f)\n", x, y);
			//cout << endl;

			intersections.push_back(pair<double, double>(x, y));
			//in_it++;
		}
		//cout << endl;
	}


	return intersections;

}

vector<double> getSphereCenterAndRadius(Image* image, const unsigned short& threshold,
                                                       const unsigned short& sphere_label) {
	vector<double> x_y_radius;

	Image shaded_image(*image);
	toBinary(&shaded_image, threshold);
	separateObjectsByShade(&shaded_image, sphere_label, 50);

	unordered_map<unsigned short, double> label_to_x_moment = getCenterX(&shaded_image);

	while(label_to_x_moment.size() > 1) {
		unsigned int t;
		shaded_image = *image;
		cout << "ERROR: your threshold value creates multiple objects.\n";
		cout << "Enter a new threshold: ";
		cin >> t;
		toBinary(&shaded_image, t);
		separateObjectsByShade(&shaded_image, sphere_label, 50);
		label_to_x_moment = getCenterX(&shaded_image);
	}

  //TODO: remove debug line
  (*image) = shaded_image;

	unordered_map<unsigned short, double> label_to_y_moment = getCenterY(&shaded_image);

	double x_center, y_center, cur_label, radius;
	for(const auto& i: label_to_x_moment) {
		//printf("(%f, %f)\n", i.second, label_to_y_moment[i.first]);
		cur_label = i.first, x_center = i.second, y_center = label_to_y_moment[i.first];
	}

	size_t radius_top = 0, radius_bottom = 0, radius_left = 0, radius_right = 0;

	size_t i = x_center, j = y_center;
	while(shaded_image.GetPixel(i, j) == sphere_label) {
		j++;
		radius_right++;
	}
	j = y_center;
	while(shaded_image.GetPixel(i, j) == sphere_label) {
		j--;
		radius_left++;
	}
  j = y_center;
  while(shaded_image.GetPixel(i, j) == sphere_label) {
    i--;
    radius_top++;
  }
  i = x_center;
  while(shaded_image.GetPixel(i, j) == sphere_label) {
    i++;
    radius_bottom++;
  }

  double average_diameter = (double)(radius_left + radius_right + radius_top + radius_bottom) / 2;
  x_y_radius.push_back(x_center);
  x_y_radius.push_back(y_center);
  x_y_radius.push_back(average_diameter/2);
  //x_y_radius_label.push_back(label);
	return x_y_radius;
}
pair<double, double> getHighlight(const Image* image, const unsigned short& threshold,
                                                const unsigned short& sphere_label) {
  //threshold to a greater value
  pair<double, double> highlight_coordinate;
  Image labeled_image(*image);

  toBinary(&labeled_image, threshold);
  unsigned short highlight_label = 255;
  separateObjectsByShade(&labeled_image, highlight_label, 0);
  unordered_map<unsigned short, double> highlight_center_x = getCenterX(&labeled_image);
  unordered_map<unsigned short, double> highlight_center_y = getCenterY(&labeled_image);

  highlight_coordinate = pair<double, double>(highlight_center_x[highlight_label],
                                              highlight_center_y[highlight_label]);
  return highlight_coordinate;

}
double getDistanceFromSphereCenter(pair<double, double> center, pair<double, double> target) {

  double radicant = (target.first-center.first)*(target.first-center.first)
                    + (target.second-center.second)*(target.second-center.second);

  return sqrt(radicant);

}
vector<double> getUnitVector(vector<double> v) {

  double length = 0;
  vector<double> unitVector;

  for(const auto& i: v) {
    length += (i*i);
  }

  length = sqrt(length);

  for(const auto& i: v) {
   unitVector.push_back(i/length);
  }
  return unitVector;

}

double getVectorLength(vector<double> v) {
  double length = 0;
  for(const auto& i: v) {
    length += (i*i);
  }
  return sqrt(length);
}

vector<vector<double>> getLightSourceMatrix(const Image* image1,
                                            const Image* image2,
                                            const Image* image3,
                                            const unsigned short& threshold1,
                                            const unsigned short& threshold2,
                                            const unsigned short& threshold3,
                                            double centerX, double centerY, double radius) {

  pair<double, double> center(centerX, centerY);
  vector<vector<double>> lightsource_matrix;

  pair<double, double> highlight1 = getHighlight(image1, threshold1, 50);
  pair<double, double> highlight2 = getHighlight(image2, threshold2, 50);
  pair<double, double> highlight3 = getHighlight(image3, threshold3, 50);
/*
  double distance1 = getDistanceFromSphereCenter(center, highlight1);
  double distance2 = getDistanceFromSphereCenter(center, highlight2);
  double distance3 = getDistanceFromSphereCenter(center, highlight3);
*/
  double x_distance1 = highlight1.first - centerX;
  double x_distance2 = highlight2.first - centerX;
  double x_distance3 = highlight3.first - centerX;

  double y_distance1 = highlight1.second - centerY;
  double y_distance2 = highlight2.second - centerY;
  double y_distance3 = highlight3.second - centerY;

  double z_distance1 = sqrt(radius*radius - x_distance1*x_distance1 - y_distance1*y_distance1);
  double z_distance2 = sqrt(radius*radius - x_distance2*x_distance2 - y_distance2*y_distance2);
  double z_distance3 = sqrt(radius*radius - x_distance3*x_distance3 - y_distance3*y_distance3);

/*
  double angle1 = acos(distance1/radius);
  double angle2 = acos(distance2/radius);
  double angle3 = acos(distance3/radius);
*/


  vector<double> lightsource_vector1 = {x_distance1,
                                        y_distance1,
                                        z_distance1};
  vector<double> lightsource_vector2 = {x_distance2,
                                        y_distance2,
                                        z_distance2};
  vector<double> lightsource_vector3 = {x_distance3,
                                        y_distance3,
                                        z_distance3};

/*
  vector<double> lightsource_vector1 = {highlight1.first-centerX,
                                        highlight1.second-centerY,
                                        radius*sin(angle1)};
  vector<double> lightsource_vector2 = {highlight2.first-centerX,
                                        highlight2.second-centerY,
                                        radius*sin(angle2)};
  vector<double> lightsource_vector3 = {highlight3.first-centerX,
                                        highlight3.second-centerY,
                                        radius*sin(angle3)};
*/
  lightsource_vector1 = getUnitVector(lightsource_vector1);
  lightsource_vector2 = getUnitVector(lightsource_vector2);
  lightsource_vector3 = getUnitVector(lightsource_vector3);

  unsigned short highlight1_value = image1->GetPixel(highlight1.first,
                                                     highlight1.second);
  unsigned short highlight2_value = image2->GetPixel(highlight2.first,
                                                     highlight2.second);
  unsigned short highlight3_value = image3->GetPixel(highlight3.first,
                                                     highlight3.second);

  for(size_t i = 0; i < lightsource_vector1.size(); ++i) {
    lightsource_vector1[i] *= highlight1_value;
  }
  for(size_t i = 0; i < lightsource_vector2.size(); ++i) {
    lightsource_vector2[i] *= highlight2_value;
  }
  for(size_t i = 0; i < lightsource_vector3.size(); ++i) {
    lightsource_vector3[i] *= highlight3_value;
  }

  lightsource_matrix.push_back(lightsource_vector1);
  lightsource_matrix.push_back(lightsource_vector2);
  lightsource_matrix.push_back(lightsource_vector3);

  return lightsource_matrix;

}

Image getNeedleImage(const Image* image1, const Image* image2, const Image* image3,
                     const unsigned short& threshold, const unsigned int& step,
                     const vector<vector<double>> lightsource_matrix) {

  Image needle_image(*image1);
  vector<double> unit_normal;
  vector<unsigned short> intensity;
  vector<vector<double>> lightsource_matrix_inverse = invert3x3Matrix(lightsource_matrix);

  for(size_t i = 0; i < needle_image.num_rows(); i += step) {
    for(size_t j = 0; j < needle_image.num_columns(); j += step) {

      if(image1->GetPixel(i, j) > threshold &&
         image2->GetPixel(i, j) > threshold &&
         image3->GetPixel(i, j) > threshold) {

         intensity = {image1->GetPixel(i, j),
                      image2->GetPixel(i, j),
                      image3->GetPixel(i, j)};

         unit_normal = calculateNormalVector(lightsource_matrix_inverse,
                                             intensity);
         unit_normal = getUnitVector(unit_normal);

         for(size_t i = 0; i < unit_normal.size(); ++i) {
           unit_normal[i] *= 10;
         }
         DrawLine(i, j, i+unit_normal[0], j+unit_normal[1], 255, &needle_image);

         needle_image.SetPixel(i, j, 0);
         if(i != 0) {
          needle_image.SetPixel(i-1, j, 255);
         }
         if(j != 0) {
          needle_image.SetPixel(i, j-1, 255);
         }
         if((i+1) < needle_image.num_rows()) {
          needle_image.SetPixel(i+1, j, 255);
         }
         if((j+1) < needle_image.num_columns()) {
          needle_image.SetPixel(i, j+1, 255);
         }

      }

    }
  }

  return needle_image;
}
Image getAlbedoImage(const Image* image1, const Image* image2, const Image* image3,
                     const unsigned short& threshold,
                     const vector<vector<double>> lightsource_matrix) {

  Image albedo_image(*image1);
  vector<double> normal;
  vector<unsigned short> intensity;
  vector<vector<double>> lightsource_matrix_inverse = invert3x3Matrix(lightsource_matrix);
  double maximum_albedo = 0.0;
  // find maximum albedo
  for(size_t i = 0; i < albedo_image.num_rows(); ++i) {
    for(size_t j = 0; j < albedo_image.num_columns(); ++j) {

      if(image1->GetPixel(i, j) > threshold &&
         image2->GetPixel(i, j) > threshold &&
         image3->GetPixel(i, j) > threshold) {

        intensity = {image1->GetPixel(i, j),
                     image2->GetPixel(i, j),
                     image3->GetPixel(i, j)};

        normal = calculateNormalVector(lightsource_matrix_inverse,
                                       intensity);
        if(maximum_albedo < getVectorLength(normal))
          maximum_albedo = getVectorLength(normal);
      }
    }
  }

  double current_albedo;

  for(size_t i = 0; i < albedo_image.num_rows(); ++i) {
    for(size_t j = 0; j < albedo_image.num_columns(); ++j) {

      if(image1->GetPixel(i, j) > threshold &&
         image2->GetPixel(i, j) > threshold &&
         image3->GetPixel(i, j) > threshold) {

        intensity = {image1->GetPixel(i, j),
                     image2->GetPixel(i, j),
                     image3->GetPixel(i, j)};

        normal = calculateNormalVector(lightsource_matrix_inverse,
                                       intensity);
        current_albedo = getVectorLength(normal);
        albedo_image.SetPixel(i, j, 255*(current_albedo/maximum_albedo));
      }
      else {
        albedo_image.SetPixel(i, j, 0);
      }
    }
  }


  return albedo_image;

}
vector<vector<double>> invert3x3Matrix(const vector<vector<double>>& matrix) {
  vector<vector<double>> inverted;

  double cofactor11 = matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1];
  double cofactor12 = -1*(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]);
  double cofactor13 = matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0];

  double determinant = matrix[0][0]*cofactor11 +
                       matrix[0][1]*cofactor12 +
                       matrix[0][2]*cofactor13;
  //cout << "determinant: " << determinant << endl;

  double cofactor21 = -1*(matrix[0][1]*matrix[2][2] - matrix[0][2]*matrix[2][1]);
  double cofactor22 = matrix[0][0]*matrix[2][2] - matrix[0][2]*matrix[2][0];
  double cofactor23 = -1*(matrix[0][0]*matrix[2][1] - matrix[0][1]*matrix[2][0]);

  double cofactor31 = matrix[0][1]*matrix[1][2] - matrix[0][2]*matrix[1][1];
  double cofactor32 = -1*(matrix[0][0]*matrix[1][2] - matrix[0][2]*matrix[1][0]);
  double cofactor33 = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];

  inverted = {{cofactor11/determinant, cofactor21/determinant, cofactor31/determinant},
              {cofactor12/determinant, cofactor22/determinant, cofactor32/determinant},
              {cofactor13/determinant, cofactor23/determinant, cofactor33/determinant}};

  return inverted;

}

vector<double> matrixMultiplySquareWithColumn(vector<vector<double>> square,
                                              vector<unsigned short> column) {
  vector<double> N;
  double temp = 0.0;

  for(const auto& i: square) {
    for(size_t j = 0; j < column.size(); ++j) {
      temp += (double)(i[j]*column[j]);
    }
    N.push_back(temp);
    temp = 0.0;
  }

  return N;
}
vector<double> calculateNormalVector(vector<vector<double>> inverse_source_direction,
                                     vector<unsigned short> intensity) {

  vector<double> albedo_normal = matrixMultiplySquareWithColumn(
                                   inverse_source_direction,
                                   intensity);

  return albedo_normal;

}

}  // namespace ComputerVisionProjects


/*

set<list<pair<size_t, size_t>>> hysteresis_threshold(vector<vector<double>>& strength,
													const vector<vector<double>>& orientation,
													double tl, double th) {


	set<list<pair<size_t, size_t>>> edge_set;
	list<pair<size_t, size_t>> edge_list;
	pair<size_t, size_t> coordinate;

	if(tl > th) {
		cout << "Invalid threshold values! Condition must be: (tl < th)" << endl;
		cout << "Returning empty set!" << endl;
		return edge_set;
	}

	vector<vector<bool>> visited(orientation.size(),
															 vector<bool>(orientation[0].size()));

	vector<vector<double>> edge_orientation(orientation.size(),
																					vector<double>(orientation[0].size(), false));
	double theta, rho, constant, tan_theta;
	long x, y;

	for(size_t i = 0; i < strength.size(); ++i) {
		for(size_t j = 0; j < strength.size(); ++j) {

			if(strength[i][j] > th && !visited[i][j]) {
				visited[i][j] = true;
				// calculate line
				// xsint - ycost + p = 0
				// the angle of the line is perpendicular to the arctan(dy, dx)
				theta = orientation[i][j];
				if(theta < 0) {
					theta = theta + M_PI;
				}
				theta = fmod_math(theta, M_PI);
				edge_orientation[i][j] = theta;
				// solve for p:
				rho = (j*cos(theta)) - (i*sin(theta));
				constant = (rho / cos(theta));
				tan_theta = tan(theta);

				x = i, y = j;
				edge_list.push_back(pair<size_t, size_t>(i, j));

				if( (theta < (M_PI/2) + 0.001) && (theta > (M_PI/2) - 0.001) ) {
					y++;
					// loop one way
					while(y < strength[0].size()) {
						if(strength[x][y] > tl) {
							coordinate.first = x;
							coordinate.second = y;
							visited[x][y] = true;
							edge_list.push_back(coordinate);
						}
						y++;
					}
					y = j-1;
					// loop the other
					while(y >= 0) {
						if(strength[x][y] > tl) {
							coordinate.first = x;
							coordinate.second = y;
							visited[x][y] = true;
							edge_list.push_back(coordinate);
						}
						y--;
					}
				} // end theta is pi/2
				else {
					x = i+1;
					y = round(x*tan_theta + constant);
					// loop one way
					while(x < strength.size() && y >= 0 && y < strength[0].size() ) {

						if(strength[x][y] > tl) {
							coordinate.first = x;
							coordinate.second = y;
							visited[x][y] = true;
							edge_list.push_back(coordinate);
						}
						x++;
						y = round(x*tan_theta + constant);
					}
					x = i-1;
					y = round(x*tan_theta + constant);
					// loop the other
					while(x >= 0 && y >= 0 && y < strength[0].size()) {

						if(strength[x][y] > tl) {
							coordinate.first = x;
							coordinate.second = y;
							visited[x][y] = true;
							edge_list.push_back(coordinate);
						}
						x--;
						y = round(x*tan_theta + constant);
					}

				}
				edge_set.insert(edge_list);
				edge_list = list<pair<size_t, size_t>>();
			}

		}

	}


}




*/
