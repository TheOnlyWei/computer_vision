SAMPLE PROGRAMS: CALCULATE ALBEDO AND NEEDLE MAP

	Type make all to make all four sample programs s1, s2, s3, s4
	Follow the DIRECTIONS to apply the programs.

		1. The formula implemented for finding normal to a sphere is (used in s2):
		(x-xc)^2 + (y-yc)^2 + (z-zc)^2 = r^2

		Where z-zc is solved to find the depth of a normal vector.

	DIRECTIONS

		s1 command line parameters:
			{input original image} {input threshold value} {output parameters file}

		s2 command line parameters:
			{input parameters file} {image 1} {image 2} {image 3} {output directions file}

		s3 command line parameters:
			{input directions} {image 1} {image 2} {image 3} {step} {threshold} {output}

		s4 command line parameters:
			{input directions} {image 1} {image 2} {image 3} {threshold} {output}

		Recommended commands with parameters:
		s1 command:
		./s1 sphere0.pgm 86 sphere_data.txt

		s2 command:
		./s2 sphere_data.txt sphere1.pgm sphere2.pgm sphere3.pgm light_direction.txt

		s3 command:
		./s3 light_direction.txt object1.pgm object2.pgm object3.pgm 10 84 needle_object.pgm
		./s3 light_direction.txt sphere1.pgm sphere2.pgm sphere3.pgm 10 84 needle_sphere.pgm

		s4 command:
		./s4 light_direction.txt object1.pgm object2.pgm object3.pgm 84 albedo_object.pgm
		./s4 light_direction.txt sphere1.pgm sphere2.pgm sphere3.pgm 84 albedo_sphere.pgm

		Copy and paste into terminal to run all programs on object files
		./s1 sphere0.pgm 86 sphere_data.txt
		./s2 sphere_data.txt sphere1.pgm sphere2.pgm sphere3.pgm light_direction.txt
		./s3 light_direction.txt object1.pgm object2.pgm object3.pgm 10 84 needle_object.pgm
		./s4 light_direction.txt object1.pgm object2.pgm object3.pgm 84 albedo_object.pgm


	Input file:
		s1: sphere0.pgm
		s2: sphere_data.txt sphere1.pgm sphere2.pgm sphere3.pgm
		s3: light_direction.txt object1.pgm object2.pgm object3.pgm
		s4: light_direction.txt object1.pgm object2.pgm object3.pgm
	Output file:
		s1: sphere_data.txt
		s2: light_direction.txt
		s3: needle_object.pgm
		s4: albedo_object.pgm
