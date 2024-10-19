/*
 * add_Input.c
 *
 *  Created on: Feb 6, 2012
 *      Author: viriyopa
 */

#include <stdlib.h>
#include <stdio.h>

#include "csln.h"


void add_Input_v2(int N_inputs_X_bef_or_eq[], int N_inputs_X_aft[], int N_inputs_EI_bef_or_eq[], int N_inputs_EI_aft[], v_pert_elem_struct *v_pert, double dvdt, double t, double N_input, input_struct *head[], input_struct *tail[], double latency, int from, int i_LE)
{
	int is_pass_t_ref_node = 0;

	/* Get memory of time variable. */
	input_struct *tmp=(input_struct *)malloc(sizeof(input_struct));
	tmp->pre		=NULL;
	tmp->arr_time	=t + latency;
	tmp->N_input	=N_input;
	tmp->from		=from;
	tmp->v_pert		=v_pert;
	tmp->dvdt		=dvdt;
	tmp->next		=NULL;

	if (tail[i_LE] == NULL)
	{
		/* The input list is empty */
		head[i_LE] = tmp;
		tail[i_LE] = tmp;
	}
	else
	{
		input_struct *head1 = head[i_LE];
		while ((head1 != NULL) && (head1->arr_time < tmp->arr_time))
		{
			if (head1->v_pert == NULL)
			{
				is_pass_t_ref_node = 1;
			}

			head1 = head1->next;
		}

		if (head1 == NULL)
		{// Reach tail of the list
			tmp->pre			=tail[i_LE];
			tail[i_LE]->next	=tmp;

			tail[i_LE]=tmp;
		}
		else
		{
			input_struct *pre_head1 = head1->pre;
			if (pre_head1 == NULL)
			{// Reach head of the list
				tmp->next = head1;
				head1->pre = tmp;

				head[i_LE]=tmp;
			}
			else
			{
				pre_head1->next = tmp;
				tmp->pre = pre_head1;

				tmp->next = head1;
				head1->pre = tmp;
			}
		}

		if (is_pass_t_ref_node == 0)
		{
			N_inputs_X_bef_or_eq[i_LE]++;
			N_inputs_EI_bef_or_eq[i_LE]++;
		}
		else
		{
			N_inputs_X_aft[i_LE]++;
			N_inputs_EI_aft[i_LE]++;
		}
	}
}

void add_Input(v_pert_elem_struct *v_pert, double dvdt, double t, double N_input, input_struct *head[], input_struct *tail[], double latency, int from, int input_i, int is_all_delays_the_same)
{
	/* Get memory of time variable. */
	input_struct *tmp=(input_struct *)malloc(sizeof(input_struct));
	tmp->pre		=NULL;
	tmp->arr_time	=t + latency;
	tmp->N_input	=N_input;
	tmp->from		=from;
	tmp->v_pert		=v_pert;
	tmp->dvdt		=dvdt;
	tmp->next		=NULL;

	if (tail[input_i] == NULL)
	{
		/* The input list is empty */
		head[input_i] = tmp;
		tail[input_i] = tmp;
	}
	else
	{
		if (is_all_delays_the_same == 0)
		{
			/* The input list is not empty */
			if (tail[input_i]->arr_time <= tmp->arr_time)
			{
				tmp->pre			=tail[input_i];
				tail[input_i]->next	=tmp;

				tail[input_i]=tmp;
			}
			else
			{
				input_struct *tmp2=NULL;
				input_struct *tail1=tail[input_i];
				while (tail1!=NULL)
				{
					if (tail1->arr_time < tmp->arr_time)
					{
						input_struct *tmp3=tail1->next;

						tmp->pre	=tail1;
						tail1->next	=tmp;

						tmp3->pre	=tmp;
						tmp->next	=tmp3;

						return;
					}

					tmp2=tail1;
					tail1=tail1->pre;
				}

				/* Append at the head of the list */
				tmp2->pre=tmp;
				tmp->next=tmp2;

				head[input_i]=tmp;
			}
		}
		else
		{
			tmp->pre			=tail[input_i];
			tail[input_i]->next	=tmp;

			tail[input_i]=tmp;
		}
	}
}


