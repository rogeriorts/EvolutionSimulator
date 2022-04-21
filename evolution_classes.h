#pragma once
#include <cmath>
#include  <iostream>




struct grid_element {
	int i;
	int j;
	const unsigned char* col;
	int* cells_id_list;
	int last_cell_id_on_grid_element;
	double energy;

};

struct recipe {
	int i_energy_sensors_positions[4];
	int j_energy_sensors_positions[4];
	double i_enegy_sensor_weights[4];
	double j_enegy_sensor_weights[4];
	double i_enegy_sensor_bias;
	double j_enegy_sensor_bias;

	int max_allowed_movement;
	double energy_absorption_efficiency;
	int max_number_of_cells;

	double reproduce_probability;
};


struct organism {

	int* cells_ids;
	int last_cell_id;

	int* organism_sensor_cells_ids;
	int last_sensor_cell_id;

	int* organism_energy_cells_ids;
	int last_energy_cell_id;

	recipe organism_recipe;

	double resultant_i_movement;
	double resultant_j_movement;
	double energy;
	bool alive;
};


struct cell {
	int i;
	int j;
	int grid_id;
	int index_of_cell_id_in_grid_list;
	int cell_type;
	int organism_id;
	bool alive;
};