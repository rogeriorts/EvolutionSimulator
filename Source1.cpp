
#include  <iostream>
#include "CImg.h"
#include <cmath>
#include "util.h"
#include <time.h>
#include "evolution_classes.h"

using namespace std;
using namespace cimg_library;




int GetGridArrayPosition(int i, int j, int number_of_squares_in_i)
{
	return j * number_of_squares_in_i + i;
}


int ComputePosition(int s, int delta_s, int number_of_squares_s)
{


	if (s + delta_s < 0)
		return number_of_squares_s + s + delta_s;
	if (s + delta_s > number_of_squares_s - 1)
		return s + delta_s - number_of_squares_s;
	else return s + delta_s;
}

void GetEnergyArray(int i, int j, recipe &organism_recipe, double* output_energy, grid_element*& grid, int number_of_squares_in_i, int number_of_squares_in_j)
{
	for (int k = 0; k < 4; k++) {

		int i_measurement = ComputePosition(i,organism_recipe.i_energy_sensors_positions[k],number_of_squares_in_i);
		int j_measurement = ComputePosition(j, organism_recipe.j_energy_sensors_positions[k], number_of_squares_in_j);

		int grid_array_position = GetGridArrayPosition(
			i_measurement,
			j_measurement,
			number_of_squares_in_i
		);

		output_energy[k] = grid[grid_array_position].energy;
		

	}
}



void CalculateResultantMovementFromSensors(organism* &organism_, int n_organisms, cell* &cells, grid_element* &grid, int number_of_squares_in_i, int number_of_squares_in_j)
{
	for (int o = 0; o < n_organisms; o++)
	{
		organism_[o].resultant_i_movement = 0.0;
		organism_[o].resultant_j_movement = 0.0;
		for (int c = 0; c <=organism_[o].last_sensor_cell_id; c++)
		{

			int i = cells[organism_[o].organism_sensor_cells_ids[c]].i;
			int j = cells[organism_[o].organism_sensor_cells_ids[c]].j;

			double energy[4];

			GetEnergyArray(i, j, organism_[o].organism_recipe, energy, grid, number_of_squares_in_i, number_of_squares_in_j);

			organism_[o].resultant_i_movement += DotProduct(energy, organism_[o].organism_recipe.i_enegy_sensor_weights, 4);
			organism_[o].resultant_j_movement += DotProduct(energy, organism_[o].organism_recipe.j_enegy_sensor_weights, 4);
			organism_[o].resultant_i_movement += organism_[o].organism_recipe.i_enegy_sensor_bias;
			organism_[o].resultant_j_movement += organism_[o].organism_recipe.j_enegy_sensor_bias;

		}
	}
}


void AddCellToGrid(grid_element * &grid,int i, int j, int number_of_squares_in_i, cell* &cells,int cell_id)
{
	

	int grid_id = number_of_squares_in_i * j + i;

	cells[cell_id].i = i;
	cells[cell_id].j = j;
	cells[cell_id].grid_id = grid_id;

	grid[grid_id].last_cell_id_on_grid_element++;

	int p = grid[grid_id].last_cell_id_on_grid_element;

	cells[cell_id].index_of_cell_id_in_grid_list = p;

	grid[grid_id].cells_id_list[p] = cell_id;

}

void RemoveCellFromGrid(grid_element* &grid, int number_of_squares_in_i, cell*& cells, int cell_id)
{
	int grid_id = cells[cell_id].grid_id;
	
	int last_cell_id_on_grid_element = grid[grid_id].last_cell_id_on_grid_element;

	int current_cell_position_in_cells_list = cells[cell_id].index_of_cell_id_in_grid_list;
	
	grid[grid_id].cells_id_list[current_cell_position_in_cells_list] = grid[grid_id].cells_id_list[last_cell_id_on_grid_element];

	grid[grid_id].cells_id_list[last_cell_id_on_grid_element] = -1;

	grid[grid_id].last_cell_id_on_grid_element--;

	cells[cell_id].i = -1;
	cells[cell_id].j = -1;
	cells[cell_id].grid_id = -1;


}


bool CanMoveInDirection(int s, int delta_s, int number_of_squares_s)
{
	if (s + delta_s < 0) return false;
	if (s + delta_s > number_of_squares_s - 1) return false;
	return true;
}




void MoveCell(grid_element* &grid, int number_of_squares_in_i, int number_of_squares_in_j, 
	int delta_i, int delta_j, cell* cells, int cell_id)
{
	int i = ComputePosition(cells[cell_id].i, delta_i, number_of_squares_in_i);
	int j = ComputePosition(cells[cell_id].j, delta_j, number_of_squares_in_j);

	RemoveCellFromGrid(grid, number_of_squares_in_i, cells, cell_id);
	AddCellToGrid(grid, i, j, number_of_squares_in_i, cells, cell_id);

}



void CreateGrid(int number_of_squares_in_i, int number_of_squares_in_j, grid_element* &grid, double energy_source_multiplicator)
{
	grid = new grid_element[number_of_squares_in_i * number_of_squares_in_j];
	
	int array_pos = 0;
	
	for (int j = 0; j < number_of_squares_in_j; j++)
		for (int i = 0; i < number_of_squares_in_i; i++)
		{
			grid[array_pos].i = i;
			grid[array_pos].j = j;
			double squared_radius = (double)(i * i + j * j);

			//grid[array_pos].energy = abs(energy_source_multiplicator * pow((((double)number_of_squares_in_i - sqrt(squared_radius))/ (double)number_of_squares_in_i),3));
			grid[array_pos].cells_id_list = new int[100];
			for (int k = 0; k < 100; k++) grid[array_pos].cells_id_list[k] = -1;
			grid[array_pos].last_cell_id_on_grid_element = -1;

			array_pos++;

		}

}

void SetEnergyToGrid(int number_of_squares_in_i, int number_of_squares_in_j, grid_element*& grid, double energy_source_multiplicator, int iter)
{
	int array_pos = 0;

	for (int j = 0; j < number_of_squares_in_j; j++)
		for (int i = 0; i < number_of_squares_in_i; i++)
		{
			grid[array_pos].i = i;
			grid[array_pos].j = j;
			double squared_radius = (double)(i * i + j * j);

			grid[array_pos].energy = abs(energy_source_multiplicator * (
					cos(2 * 3.1416 / number_of_squares_in_i *  (i+ 100*cos((double)iter/50.0))
				)
			));
			array_pos++;

		}

}

void DestroyGrid(int number_of_squares_in_i, int number_of_squares_in_j, grid_element*& grid)
{
	

	for (int array_pos = 0; array_pos < number_of_squares_in_i* number_of_squares_in_j; array_pos++)		
	{
		delete[] grid[array_pos].cells_id_list;
	}
	delete[] grid;
}


void MoveOrganism(cell* &cells, organism* &organisms, int organism_id, int delta_i,int delta_j, int number_of_squares_in_i, int number_of_squares_in_j,  grid_element*& grid)
{

	//eventually, check if things are mobile
	for ( int c = 0; c <= organisms[organism_id].last_cell_id; c++)
	{
		int cell_id = organisms[organism_id].cells_ids[c];
		MoveCell(grid, number_of_squares_in_i, number_of_squares_in_j, delta_i, delta_j, cells, cell_id);
	}
		
}

void OrganismsGetEnergyFromEnvironmnet(cell*& cells, organism*& organisms, int n_organisms, int number_of_squares_in_i, grid_element*& grid)
{
	for (int organism_id = 0; organism_id < n_organisms; organism_id++)
		if(organisms[organism_id].alive)
		{
			for (int c = 0; c <= organisms[organism_id].last_energy_cell_id; c++)
			{
				int cell_id = organisms[organism_id].organism_energy_cells_ids[c];
				int grid_id = GetGridArrayPosition(cells[cell_id].i, cells[cell_id].j, number_of_squares_in_i);
				double env_energy = grid[grid_id].energy;



				organisms[organism_id].energy += (grid[grid_id].energy * organisms[organism_id].organism_recipe.energy_absorption_efficiency);

				if (organisms[organism_id].energy > 5 * organisms[organism_id].last_cell_id + 1)
					organisms[organism_id].energy = 5 * organisms[organism_id].last_cell_id + 1;

			}
		}
}


void MoveOrganisms(cell*& cells, organism*& organisms, int n_organisms, int number_of_squares_in_i, int number_of_squares_in_j, grid_element*& grid)
{

	
	CalculateResultantMovementFromSensors(organisms, n_organisms, cells, grid, number_of_squares_in_i, number_of_squares_in_j);

	for (int organism_id=0;organism_id<n_organisms;organism_id++)
	{
		if(organisms[organism_id].alive)
		{
			int sgn_i = sgn(organisms[organism_id].resultant_i_movement);
			int sgn_j = sgn(organisms[organism_id].resultant_j_movement);

			double abs_i = abs(organisms[organism_id].resultant_i_movement);
			double abs_j = abs(organisms[organism_id].resultant_j_movement);

			int delta_i = sgn_i * min((int) round(abs_i), organisms[organism_id].organism_recipe.max_allowed_movement);
			int delta_j = sgn_j * min((int) round(abs_j), organisms[organism_id].organism_recipe.max_allowed_movement);
		
			int n_cells_in_organism = 0;

			for (int c = 0; c <= organisms[organism_id].last_cell_id; c++)
			{
				int cell_id = organisms[organism_id].cells_ids[c];
				MoveCell(grid, number_of_squares_in_i, number_of_squares_in_j, delta_i, delta_j, cells, cell_id);
			}
			organisms[organism_id].energy -= organisms[organism_id].last_energy_cell_id + organisms[organism_id].last_sensor_cell_id + 2+abs_i+abs_j;

		}
	}
}


void OrganismsDieFromStarvation(cell*& cells, organism*& organisms, int n_organisms, int & alive_organisms)
{
	for (int organism_id = 0; organism_id < n_organisms; organism_id++)
	{
		if (organisms[organism_id].alive)
		{
			if (organisms[organism_id].energy <= 0)
			{
				organisms[organism_id].alive = false;
				alive_organisms--;
				for ( int c = 0; c <= organisms[organism_id].last_cell_id; c++)
				{
					int cell_id = organisms[organism_id].cells_ids[c];
					cells[cell_id].alive = false;


				}

			}
		}


	}
}







void CreateOrganism(cell*& cells, organism*& organisms, int n_cells, int n_cells_max_per_organism, int i_max, int j_max, int i_min, int j_min, int number_of_squares_in_i, grid_element*& grid, int& last_org_id, int n_org_max)
{

	if (last_org_id < n_org_max)
	{
		last_org_id++;
		int o = last_org_id;
		//create organism recipe
		organisms[o].organism_recipe.max_number_of_cells = n_cells_max_per_organism;
		organisms[o].organism_recipe.max_allowed_movement = GenerateIntegerRandomNumber(1, 2);
		organisms[o].organism_recipe.energy_absorption_efficiency = GenerateRandomNumber(0.0, 1.0);
		organisms[o].organism_recipe.reproduce_probability = GenerateRandomNumber(0.0, 1.0);

		for (int r = 0; r < 4; r++)
		{
			organisms[o].organism_recipe.i_enegy_sensor_weights[r] = GenerateRandomNumber(-1.0, 1.0);
			organisms[o].organism_recipe.j_enegy_sensor_weights[r] = GenerateRandomNumber(-1.0, 1.0);
			organisms[o].organism_recipe.i_energy_sensors_positions[r] = GenerateIntegerRandomNumber(-50, 50);
			organisms[o].organism_recipe.j_energy_sensors_positions[r] = GenerateIntegerRandomNumber(-50, 50);
		}

		organisms[o].organism_recipe.i_enegy_sensor_bias = GenerateRandomNumber(-1.0, 1.0);
		organisms[o].organism_recipe.j_enegy_sensor_bias = GenerateRandomNumber(-1.0, 1.0);

		//initializing arrays
		organisms[o].last_cell_id = -1;
		organisms[o].last_sensor_cell_id = -1;
		organisms[o].last_energy_cell_id = -1;

		//create_cells

		organisms[o].cells_ids = new int[n_cells_max_per_organism];
		organisms[o].organism_sensor_cells_ids = new int[n_cells_max_per_organism];
		organisms[o].organism_energy_cells_ids = new int[n_cells_max_per_organism];


		//placing organisms_cells
		int i = GenerateIntegerRandomNumber(i_min, i_max);
		int j = GenerateIntegerRandomNumber(j_min, j_max);

		int i_growth = GenerateRandomNumber(0.0, 1.0) >= 0.5 ? 1 : -1;
		int j_growth = GenerateRandomNumber(0.0, 1.0) >= 0.5 ? 1 : -1;

		for (int c = 0; c < n_cells; c++)
		{
			int new_org_cell_id = o * n_cells_max_per_organism + c;
			

			organisms[o].cells_ids[c] = new_org_cell_id;
			organisms[o].last_cell_id = c;//I've got to change the name of these last_* variables, its getting confusing

			AddCellToGrid(grid, i, j, number_of_squares_in_i, cells, new_org_cell_id);

			cells[new_org_cell_id].organism_id = o;
			cells[new_org_cell_id].cell_type = (GenerateRandomNumber(0.0, 1.0) >= 0.5) ? 0 : 1;
			cells[new_org_cell_id].alive = true;

			if (cells[new_org_cell_id].cell_type == 0)//energy_cells;
			{
				organisms[o].last_energy_cell_id++;
				int index = organisms[o].last_energy_cell_id;
				organisms[o].organism_energy_cells_ids[index] = new_org_cell_id;

			}
			else if (cells[new_org_cell_id].cell_type == 1)//energy_cells;
			{
				organisms[o].last_sensor_cell_id++;
				int index = organisms[o].last_sensor_cell_id;
				organisms[o].organism_sensor_cells_ids[index] = new_org_cell_id;
			}


			if (GenerateRandomNumber(0.0, 1.0) >= 0.5)
			{
				i += i_growth;
				j_growth *= -1;
			}

			j += j_growth;

		}

		organisms[o].energy = GenerateRandomNumber(3.0, 5.0);
		organisms[o].alive = true;
		
	}
}


void CreateOrganisms(cell*& cells, organism*& organisms, int n_organisms, int n_cells, int n_cells_max_per_organism, int i_max, int j_max, int i_min, int j_min, int number_of_squares_in_i, grid_element*& grid, int& last_org_id, int n_org_max)
{


	while(true)
	{


		CreateOrganism(cells, organisms, n_cells, n_cells_max_per_organism, i_max, j_max, i_min, j_min, number_of_squares_in_i, grid, last_org_id, n_org_max);
		if (last_org_id >= n_organisms - 1) break;


	}

}




double GetClonedDouble(double src_var, double mutation_prob, double range_inf, double range_sup)
{
	if (GenerateRandomNumber(0.0, 1.0) <= mutation_prob)
		return GenerateRandomNumber(range_inf, range_sup);
	
	return src_var;
}

int GetClonedInteger(int src_var, double mutation_prob, int range_inf, int range_sup)
{
	if (GenerateRandomNumber(0.0, 1.0) <= mutation_prob)
		return GenerateIntegerRandomNumber(range_inf, range_sup);

	return src_var;
}


void CloneOrganism(cell*& cells, organism*& organisms, int source_org, int number_of_squares_in_i, int number_of_squares_in_j, grid_element*& grid,  int& last_org_id, int n_org_max, double mutation_prob,int& alive_organisms)
{

	if (organisms[source_org].energy > 2 * (organisms[source_org].last_cell_id+1))
	{
		int o = source_org;
		int new_org_id = -1;
		
		bool create_new_dynamic_arrays = false;
		
		if (last_org_id < n_org_max - 1)
		{
			last_org_id++;
			new_org_id = last_org_id;
			create_new_dynamic_arrays = true;
		}
		else
		{
			bool can_clone = false;
			for (int k = 0; k < n_org_max; k++)
			{
				if (!organisms[k].alive)
				{
					can_clone = true;
					new_org_id = k;
					break;
				}
			}
			if (!can_clone) return;
		}
		
		

		alive_organisms++;
		
		//create organism recipe
		organisms[new_org_id].organism_recipe.max_number_of_cells = organisms[o].organism_recipe.max_number_of_cells;
		
		organisms[new_org_id].organism_recipe.max_allowed_movement = GetClonedInteger(
			organisms[o].organism_recipe.max_allowed_movement,mutation_prob,1, 2);
		
		organisms[new_org_id].organism_recipe.energy_absorption_efficiency = GetClonedDouble(
			organisms[o].organism_recipe.energy_absorption_efficiency, mutation_prob, 0.0, 1.0);

		organisms[new_org_id].organism_recipe.reproduce_probability = GetClonedDouble(
			organisms[o].organism_recipe.reproduce_probability, mutation_prob, 0.0, 1.0);

		for (int r = 0; r < 4; r++)
		{
			organisms[new_org_id].organism_recipe.i_enegy_sensor_weights[r] = GetClonedDouble(
				organisms[o].organism_recipe.i_enegy_sensor_weights[r], mutation_prob, -1.0, 1.0);

			organisms[new_org_id].organism_recipe.j_enegy_sensor_weights[r] = GetClonedDouble(
				organisms[o].organism_recipe.j_enegy_sensor_weights[r], mutation_prob, -1.0, 1.0);


			organisms[new_org_id].organism_recipe.i_energy_sensors_positions[r] = GetClonedInteger(
				organisms[o].organism_recipe.i_energy_sensors_positions[r], mutation_prob, -50, 50);


			organisms[new_org_id].organism_recipe.j_energy_sensors_positions[r] = GetClonedInteger(
				organisms[o].organism_recipe.j_energy_sensors_positions[r], mutation_prob, -50, 50);

		}


		organisms[new_org_id].organism_recipe.i_enegy_sensor_bias = GetClonedDouble(
			organisms[o].organism_recipe.i_enegy_sensor_bias, mutation_prob, -1.0, 1.0);

		organisms[new_org_id].organism_recipe.j_enegy_sensor_bias = GetClonedDouble(
			organisms[o].organism_recipe.j_enegy_sensor_bias, mutation_prob, -1.0, 1.0);
		//copy cells
		organisms[new_org_id].last_cell_id = -1;
		organisms[new_org_id].last_sensor_cell_id = -1;
		organisms[new_org_id].last_energy_cell_id = -1;

		int n_cells_max_per_organism = organisms[new_org_id].organism_recipe.max_number_of_cells;

		if(create_new_dynamic_arrays)
		{
			organisms[new_org_id].cells_ids = new int[n_cells_max_per_organism];
			organisms[new_org_id].organism_sensor_cells_ids = new int[n_cells_max_per_organism];
			organisms[new_org_id].organism_energy_cells_ids = new int[n_cells_max_per_organism];
		}

		int n_cells = organisms[o].last_cell_id + 1;

		int offset_i = (int) round(GenerateRandomNumber(0.0,1.0) * (double)n_cells);
		int offset_j = (int) round(GenerateRandomNumber(0.0, 1.0) * (double)n_cells);

		for (int c = 0; c < n_cells; c++)
		{
			int new_org_cell_id = new_org_id * n_cells_max_per_organism + c;
			int source_cell_id = source_org * n_cells_max_per_organism + c;

			organisms[new_org_id].cells_ids[c] = new_org_id * n_cells_max_per_organism + c;
			organisms[new_org_id].last_cell_id = c;//I've got to change the name of these last_* variables, its getting confusing

			//copy position and type, which is basically the only info we got to clone
			
			

			

			int i = ComputePosition(cells[source_cell_id].i, offset_i,number_of_squares_in_i);
			int j = ComputePosition(cells[source_cell_id].j, offset_j, number_of_squares_in_j);
			cells[new_org_cell_id].cell_type = cells[source_cell_id].cell_type;

			
			AddCellToGrid(grid, i, j, number_of_squares_in_i, cells, new_org_cell_id);

			cells[new_org_cell_id].organism_id = new_org_id;
			
			cells[new_org_cell_id].alive = true;

			if (cells[new_org_cell_id].cell_type == 0)//energy_cells;
			{
				organisms[new_org_id].last_energy_cell_id++;
				int index = organisms[new_org_id].last_energy_cell_id;
				organisms[new_org_id].organism_energy_cells_ids[index] = new_org_cell_id;

			}
			else if (cells[new_org_cell_id].cell_type == 1)//energy_cells;
			{
				organisms[new_org_id].last_sensor_cell_id++;
				int index = organisms[new_org_id].last_sensor_cell_id;
				organisms[new_org_id].organism_sensor_cells_ids[index] = new_org_cell_id;
			}

		}

		organisms[new_org_id].energy = organisms[source_org].energy / 2;
		organisms[source_org].energy = organisms[new_org_id].energy;
		organisms[new_org_id].alive = true;
	}
}

void ReproduceOrganisms(cell*& cells, organism*& organisms, int n_organisms, int number_of_squares_in_i, int number_of_squares_in_j, grid_element*& grid, int& last_org_id, int n_org_max, double mutation_prob, int&alive_organisms)
{

	int cell_id = -1;
	for (int o = 0; o < n_organisms; o++)
	{


		CloneOrganism(cells, organisms, o, number_of_squares_in_i, number_of_squares_in_j, grid, last_org_id, n_org_max, mutation_prob, alive_organisms);



	}

}

void DestroyOrganism(organism*& organisms, int org_id)
{


	delete[] organisms[org_id].cells_ids;
	delete[] organisms[org_id].organism_sensor_cells_ids;
	delete[] organisms[org_id].organism_energy_cells_ids;

}

void DestroyAllOrganisms(organism*& organisms, int n_organisms)
{
	for (int o = 0; o < n_organisms; o++)
	{

		DestroyOrganism(organisms, o);
	}
}



class Environment {

};


int main(int argc, const char** argv)
{
	cout << "Hello world!" << endl;
	
	unsigned int square_size = 2;
	unsigned int number_of_squares_in_i = 500;
	unsigned int number_of_squares_in_j = 500;

	unsigned int w = square_size * number_of_squares_in_i;
	unsigned int h = square_size * number_of_squares_in_j;
		 
	const unsigned char bluegreen[] = {0,170,255};
	const unsigned char orange[] = { 255,140,0 };
	const unsigned char white[] = { 255,255,255 };
	const unsigned char black[] = { 0,0,0 };
	const unsigned char red[] = { 255,0,0 };

	CImg<unsigned char> bg(w,h,1,3,255);

	bg.draw_rectangle(0, 0,w,h, bluegreen);

	CImgDisplay dsp(w, h, "Hello", 0);
	
	dsp.display(bg);

	int iter = 0;
	CImg<unsigned char> img(bg);

	grid_element *grid = nullptr;
	srand(time(NULL)); // randomize seed


	int n_cells_per_organism = 7;
	int n_organisms = 900;
	int alive_organisms = n_organisms;

	int n_organisms_max = 1000;
	
	int n_cells_max_per_organism = n_cells_per_organism;
	int n_cells = n_cells_max_per_organism * n_organisms_max;
	
	organism * organisms  = new organism[n_organisms_max];
	cell* cells = new cell[n_cells];

	int cell_id = 0;
	int last_cell_id = -1;
	int last_org_id = -1;
	double mutation_prob = 0.05;
	double energy_source_multiplicator = 3;
	
	CreateGrid(number_of_squares_in_i, number_of_squares_in_j, grid, energy_source_multiplicator);
	SetEnergyToGrid(number_of_squares_in_i, number_of_squares_in_j, grid, energy_source_multiplicator, iter);
	CreateOrganisms(cells, organisms, n_organisms, n_cells_per_organism, n_cells_max_per_organism, number_of_squares_in_i - 3,
		number_of_squares_in_j - 3, 3, 3, number_of_squares_in_i,
		grid, last_org_id,n_organisms_max);
	

	while (!dsp.is_closed() &&!dsp.is_keyESC())
	{
		img = bg;
		
		
		

		int delta_i = 0;
		int delta_j = 0;

		if (dsp.is_keyARROWDOWN())
			delta_j += 1;

		if (dsp.is_keyARROWRIGHT())
			delta_i += 1;

		if (dsp.is_keyARROWUP())
			delta_j -= 1;

		if (dsp.is_keyARROWLEFT())
			delta_i -= 1;

		//

		//if(dsp.is_key())
		if(1==1)
		{
		
			

			n_organisms = last_org_id + 1;
			//MoveOrganism(cells, organisms, 0, delta_i, delta_j, number_of_squares_in_i, number_of_squares_in_j, grid);
		
			OrganismsGetEnergyFromEnvironmnet(cells, organisms, n_organisms, number_of_squares_in_i, grid);
			MoveOrganisms(cells, organisms, n_organisms, number_of_squares_in_i, number_of_squares_in_j, grid);
			OrganismsDieFromStarvation(cells, organisms, n_organisms, alive_organisms);
			
			ReproduceOrganisms(cells, organisms, n_organisms, number_of_squares_in_i, number_of_squares_in_j, grid, last_org_id, n_organisms_max,mutation_prob,alive_organisms);

			for(int k = 0; k<number_of_squares_in_i*number_of_squares_in_j;k++)
			{
				int x = grid[k].i * square_size;
				int y = grid[k].j * square_size;

				unsigned char color[3];


				for (int channel = 0; channel < 3;channel++) 
					color[channel] = (unsigned char)((grid[k].energy)/(energy_source_multiplicator) * 255.0);

				img.draw_rectangle(x, y, x + square_size, y + square_size, color);
	
			}
			cout << alive_organisms << endl;

			for (int k = 0; k < n_cells; k++)
			{
				if (cells[k].alive)
				{
					int x = cells[k].i * square_size;
					int y = cells[k].j * square_size;

					const unsigned char *col = orange;

					if (cells[k].cell_type == 0)
						col = black;
					if (cells[k].cell_type == 1)
						col = orange;

					img.draw_rectangle(x, y, x + square_size, y + square_size, col);
			
					if (cells[k].cell_type == 0)
						for (int r = 0; r < 4; r++) {
				
							int x = (cells[k].i + organisms[cells[k].organism_id].organism_recipe.i_energy_sensors_positions[r])*square_size;
							int y = (cells[k].j + organisms[cells[k].organism_id].organism_recipe.j_energy_sensors_positions[r])* square_size;
				
							//img.draw_rectangle(x, y, x + square_size, y + square_size, red);
				

						}
				}
			}
		
		dsp.display(img);
		}

		//dsp.wait();
		iter++;
		SetEnergyToGrid(number_of_squares_in_i, number_of_squares_in_j, grid, energy_source_multiplicator, iter);
		Sleep(50);
		if (alive_organisms <= 0) break;
	}


	DestroyGrid(number_of_squares_in_i, number_of_squares_in_j, grid);
	DestroyAllOrganisms(organisms, last_org_id);
	
	delete[] cells;
	delete[] organisms;

	return 0;
}
