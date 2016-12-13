// by Wei Shi
// last modified 11/10/16
// s3

#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;

vector<vector<double>> readFromFile(const string& file_name);
void writeToFile(const string& file_name,
                 const vector<vector<double>>& direction_data);


int main(int argc, char **argv){

  if (argc!=8) {
    printf("Usage: %s <input directions> <image1> <image2> <image3> <step> <threshold> <output>\n", argv[0]);
    return 0;
  }

  const string input_directions(argv[1]);
  const string image1_file(argv[2]);
  const string image2_file(argv[3]);
  const string image3_file(argv[4]);
  const string step_str(argv[5]);
  const string threshold_str(argv[6]);
  const string output(argv[7]);

  const unsigned int step = stoi(step_str);
  const unsigned int threshold = stoi(threshold_str);

  Image image1, image2, image3;
  Image needle_image;
  if (!ReadImage(image1_file, &image1)) {
    cout <<"Can't open file " <<image1_file << endl;
    return 0;
  }
  if (!ReadImage(image2_file, &image2)) {
    cout <<"Can't open file " << image2_file << endl;
    return 0;
  }
  if (!ReadImage(image3_file, &image3)) {
    cout <<"Can't open file " << image3_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program outputs a needle image of " << output << ". The light\n";
  cout << "source vector is not normalized and is mapped into range [0,255].\n";
  cout << "Recommended pixel step interval: 10.\n";
  cout << "Recommended pixel threshold: 84.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";

  vector<vector<double>> lightsource_matrix = readFromFile(input_directions);

  needle_image = getNeedleImage(&image1, &image2, &image3, threshold, step,
                                lightsource_matrix);

  if (!WriteImage(output, needle_image)){
    cout << "Can't write to file " << output << endl;
    return 0;
  }


	return 0;
}
vector<vector<double>> readFromFile(const string& file_name) {
  vector<vector<double>> lightsource_matrix;
  fstream in;
  in.open(file_name, fstream::in);

  double input;
  vector<double> temp(3, 0.0);

  stringstream ss;
  string current_row;

  while(getline(in, current_row)) {
    ss << current_row;
    ss >> temp[0] >> temp[1] >> temp[2];
    lightsource_matrix.push_back(temp);
    temp = vector<double>(3, 0.0);
  }
  return lightsource_matrix;

}
