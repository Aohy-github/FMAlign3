// This is a demo program for showing how to call SACA_K.

#include <cstring>
#include <biovoltron/algo/sort/kiss1_sorter.hpp>
#include <span>
#include<stdio.h>
#include<stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "bwt.h"
#include<stdint.h>
#include<ctype.h>
#include <chrono>
#include <thread>
///#include <nmmintrin.h>

bwt_index bitmapper_index_params;

unsigned int debug_2 = 0;

unsigned int debug_total = 0;






void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n,
           unsigned int K, unsigned int m, int level);

inline void bwt_direct_get_sa_interval_long_back_up_1
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int j_sp, j_ep;
	unsigned int l;
	unsigned int i = 0;

	unsigned int last_sp, last_ep;

	unsigned int actually_line_sp, actually_line_ep;

	unsigned int ans_sp, ans_ep, flag;


	SA_flag_string_type tmp_SA_pop_count;

	flag = 0;


	l = sp;
	///���ǳ���224������256������ٳ���64�������ȷ�������ĸ�block
	actually_line_sp = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	///���Ƕ�224ȡ�࣬�����ȷ��block�ڲ�����Ҫɨ���λ��
	last_sp = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;



	l = ep;
	///���ǳ���224������256������ٳ���64�������ȷ�������ĸ�block
	actually_line_ep = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	///���ǳ���224������256������ٳ���64�������ȷ�������ĸ�block
	last_ep = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;




	///��˵��sp��ep��ͬһ��block�ڲ�,��count�������Բ����ص�
	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}








	ans_sp = SA_flag[actually_line_sp] >> bitmapper_index_params.SA_counter_shift_length;

	if (last_sp != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length��, ��������j=0, �����Ӧ�ð�j = SA_counter_length
		j_sp = 0;

		tmp_SA_pop_count = SA_flag[actually_line_sp] & bitmapper_index_params.SA_pop_count_mode;

		while (j_sp + bitmapper_index_params.SA_flag_warp_number <= last_sp)
		{

			ans_sp = ans_sp + __builtin_popcountll(tmp_SA_pop_count);

			j_sp = j_sp + bitmapper_index_params.SA_flag_warp_number;

			actually_line_sp++;
			tmp_SA_pop_count = SA_flag[actually_line_sp];

		}

		///����ж�����ò����sp��ep������ͬһ��block
		///����ep��ֵò�ƿ���ֱ�Ӽ̳�
		if (flag == 1)
		{
			j_ep = j_sp;

			ans_ep = ans_sp;

			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		last_sp = last_sp - j_sp;


		if (last_sp != 0)
		{
			ans_sp = ans_sp +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_sp))));
		}


	}




	(*start) = ans_sp;





	if (last_ep != 0)
	{

		if (flag==2)
		{

			tmp_SA_pop_count = SA_flag[actually_line_ep];

			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}
		}
		else
		{


			ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;

			///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
			j_ep = 0;

			tmp_SA_pop_count = SA_flag[actually_line_ep] & bitmapper_index_params.SA_pop_count_mode;


			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}

		}







	}
	else
	{
		ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;
	}






	(*length) = ans_ep - (*start);

}


inline void bwt_direct_get_sa_interval_long
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int last_sp, last_ep;

	unsigned int actually_line_sp, actually_line_ep;

	unsigned int ans_sp, ans_ep;

	unsigned int j_sp, j_ep;

	unsigned int flag = 0;

	SA_flag_string_type tmp_SA_pop_count;


	actually_line_sp = ((sp / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	last_sp = sp % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	actually_line_ep = ((ep / bitmapper_index_params.compress_SA_flag) << 8) >> 6;
	last_ep = ep % bitmapper_index_params.compress_SA_flag + SA_counter_length;



	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}

	///��actually_line_sp == actually_line_ep
	///�п���last_sp=0�� ��last_ep != 0
	///��������last_sp!=00�� ��last_ep == 0
	ans_sp = SA_flag[actually_line_sp] >> bitmapper_index_params.SA_counter_shift_length;






	if (last_sp != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		j_sp = 0;

		tmp_SA_pop_count = SA_flag[actually_line_sp] & bitmapper_index_params.SA_pop_count_mode;




		while (j_sp + bitmapper_index_params.SA_flag_warp_number <= last_sp)
		{

			ans_sp = ans_sp + __builtin_popcountll(tmp_SA_pop_count);

			j_sp = j_sp + bitmapper_index_params.SA_flag_warp_number;

			actually_line_sp++;
			tmp_SA_pop_count = SA_flag[actually_line_sp];

		}


		if (flag == 1)
		{
			ans_ep = ans_sp;

			j_ep = j_sp;

			actually_line_ep = actually_line_sp;

			flag = 2;
		}



		last_sp = last_sp - j_sp;


		if (last_sp != 0)
		{
			ans_sp = ans_sp +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_sp))));
		}


	}




	(*start) = ans_sp;










	if (last_ep != 0)
	{

		if (flag == 2)
		{

			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}



		}
		else
		{
			ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;


			///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
			j_ep = 0;
			tmp_SA_pop_count = SA_flag[actually_line_ep] & bitmapper_index_params.SA_pop_count_mode;




			while (j_ep + bitmapper_index_params.SA_flag_warp_number <= last_ep)
			{

				ans_ep = ans_ep + __builtin_popcountll(tmp_SA_pop_count);

				j_ep = j_ep + bitmapper_index_params.SA_flag_warp_number;

				actually_line_ep++;
				tmp_SA_pop_count = SA_flag[actually_line_ep];

			}


			last_ep = last_ep - j_ep;


			if (last_ep != 0)
			{
				ans_ep = ans_ep +
					__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
					& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last_ep))));
			}

		}






	}
	else
	{
		ans_ep = SA_flag[actually_line_ep] >> bitmapper_index_params.SA_counter_shift_length;
	}



	if(ans_ep > *start)
	{
		(*length) = ans_ep - (*start);
	}
	else
	{
		(*length) = 0;
	}
}




inline void bwt_direct_get_sa_interval_long_back_up
(SA_flag_string_type* SA_flag, unsigned int sp, unsigned int ep, unsigned int* start, unsigned int* length)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	SA_flag_string_type tmp_SA_pop_count;


	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;


	///����SAҪ��
	///last = l % compress_SA_flag;
	///last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///����SAҪ��
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;



	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}




	(*start) = ans;

	l = ep;



	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/
	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;



	///����SAҪ��
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///����SAҪ��
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;

	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}




	(*length) = ans - (*start);

}


inline unsigned int bwt_get_sa_restrict_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	///if (line == bitmapper_index_params.shapline) return 0;
	if (line == bitmapper_index_params.shapline)
	{
		(*accessed_sa) = 0;
		return 1;
	}



	unsigned int actually_line;
	unsigned int last;

	/**
	///����SAҪ��
	unsigned int actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;


	///����SAҪ��
	///unsigned int last = l % compress_SA_flag;
	unsigned int last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;


	actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	actually_line = actually_line + (last >> 6);




	/**
	while (((SA_flag[actually_line] << (last % bitmapper_index_params.SA_flag_warp_number))
	&bitmapper_index_params.mode_SA_flag) == 0)
	**/
	while (((SA_flag[actually_line] << (last & 63))&bitmapper_index_params.mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}



		/**
		///�ƺ�ֻ������Ҫ�ģ���Ϊ�����漰��ȡBWT��
		///���������޸�
		///actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;
		actually_line = ((l / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) + l % bitmapper_index_params.compress_occ + 32;
		delta = (bwt[actually_line / bitmapper_index_params.bwt_warp_number]
		>> ((bitmapper_index_params.bwt_warp_number - actually_line%bitmapper_index_params.bwt_warp_number - 1) * 2))
		&(bwt_string_type)3;
		**/

		actually_line = ((l >> 8)*bitmapper_index_params.acctuall_bwt_gap) + (l & 255) + 32;

		/**
		delta = (bwt[(actually_line >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) * 2))
			&(bwt_string_type)3;
		**/

		delta = (bwt[(actually_line >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) << 1))
			&(bwt_string_type)3;







		l = find_occ_fm_index(l, delta, bwt, high_occ_table);


		i++;

		if (l == bitmapper_index_params.shapline)
		{
			(*accessed_sa) = i;

			return 1;
		}






		/**
		///����SAҪ��
		///actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
		actually_line = ((l / bitmapper_index_params.compress_SA_flag)*bitmapper_index_params.acctuall_SA_flag_gap)
		/ bitmapper_index_params.SA_flag_warp_number;

		///����SAҪ��
		///last = l % compress_SA_flag;
		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;

		**/
		actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + (last >> 6);

	}

	/**
	actually_line
	= ((l / bitmapper_index_params.compress_SA_flag)*
	bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/

	actually_line = actually_line - (last >> 6);


	///����SAҪ��
	///unsigned int ans = SA_flag[actually_line];
	unsigned int ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;


	///����SAҪ��
	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;
		SA_flag_string_type tmp_SA_pop_count;
		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}

	(*accessed_sa) = (sa[ans] & bitmapper_index_params.SA_header_mode)*bitmapper_index_params.compress_sa + i;

	return 1;
}



inline unsigned int bwt_get_sa_restrict_steps_less_than_4
(SA_flag_string_type* SA_flag, unsigned int line,
bwt_string_type *bwt, high_occ_table_type* high_occ_table,
unsigned int *sa, unsigned int need_step, unsigned int* accessed_sa)
{
	unsigned int l = line;
	bwt_string_type delta;
	unsigned int i = 0;


	if (line == bitmapper_index_params.shapline) return 0;
	unsigned int actually_line;
	unsigned int last;

	/**
	///����SAҪ��
	unsigned int actually_line = ((l / bitmapper_index_params.compress_SA_flag)
	*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;


	///����SAҪ��
	///unsigned int last = l % compress_SA_flag;
	unsigned int last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;


	actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	actually_line = actually_line + (last >> 6);




	/**
	while (((SA_flag[actually_line] << (last % bitmapper_index_params.SA_flag_warp_number))
	&bitmapper_index_params.mode_SA_flag) == 0)
	**/
	while (((SA_flag[actually_line] << (last & 63))&bitmapper_index_params.mode_SA_flag) == 0)
	{

		if (i >= need_step)
		{
			return 0;
		}



		/**
		///�ƺ�ֻ������Ҫ�ģ���Ϊ�����漰��ȡBWT��
		///���������޸�
		///actually_line = ((l / compress_occ)*acctuall_bwt_gap) + l % compress_occ + 64;
		actually_line = ((l / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) + l % bitmapper_index_params.compress_occ + 32;
		delta = (bwt[actually_line / bitmapper_index_params.bwt_warp_number]
		>> ((bitmapper_index_params.bwt_warp_number - actually_line%bitmapper_index_params.bwt_warp_number - 1) * 2))
		&(bwt_string_type)3;
		**/

		actually_line = ((l >> 8)*bitmapper_index_params.acctuall_bwt_gap) + (l & 255) + 32;

		/**
		delta = (bwt[(actually_line >> 5)]
		>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) * 2))
		&(bwt_string_type)3;
		**/

		delta = (bwt[(actually_line >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line & 31) - 1) << 1))
			&(bwt_string_type)3;







		l = find_occ_fm_index(l, delta, bwt, high_occ_table);


		i++;

		if (l == bitmapper_index_params.shapline)
		{
			(*accessed_sa) = i;

			return i;
		}




		/**
		///����SAҪ��
		///actually_line = ((l / compress_SA_flag)*acctuall_SA_flag_gap) / SA_flag_warp_number + 1;
		actually_line = ((l / bitmapper_index_params.compress_SA_flag)*bitmapper_index_params.acctuall_SA_flag_gap)
		/ bitmapper_index_params.SA_flag_warp_number;

		///����SAҪ��
		///last = l % compress_SA_flag;
		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + last / bitmapper_index_params.SA_flag_warp_number;

		**/
		actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;

		last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

		actually_line = actually_line + (last >> 6);

	}

	/**
	actually_line
	= ((l / bitmapper_index_params.compress_SA_flag)*
	bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/

	actually_line = actually_line - (last >> 6);


	///����SAҪ��
	///unsigned int ans = SA_flag[actually_line];
	unsigned int ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;


	///����SAҪ��
	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;
		SA_flag_string_type tmp_SA_pop_count;
		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}

	(*accessed_sa) = sa[ans] + i;

	return 1;
}






inline void bwt_accesss_SA_cur_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step)
{

	unsigned int t;
	for (t = sp; t<ep; t++)
	{
		///����SAҪ��
		if (bwt_get_sa_restrict_steps_more_than_3(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;
		}

	}

}




inline void bwt_accesss_SA_cur_less_than_4(unsigned int *sa, SA_flag_string_type* SA_flag, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length, unsigned int need_step)
{

	unsigned int t;
	for (t = sp; t<ep; t++)
	{
		///����SAҪ��
		if (bwt_get_sa_restrict_steps_less_than_4(SA_flag, t, bwt, high_occ_table, sa,
			need_step, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;
		}

	}

}








void indenpendent_get_sa_fromFILE(unsigned int **sa, unsigned int cc, char *refer)
{
	unsigned int i;
	unsigned int n;
	n = cc-1;

	fprintf(stdout, "r_n=%u\n", n);

	n++; // append the virtual sentinel
	fprintf(stdout, "Allocating input and output space: %u bytes = %.2lf MB", 5 * n, (double)5 * n / 1024 / 1024);
	unsigned char *s_ch = new unsigned char[n];
	//unsigned char *s_ch;
	unsigned int *SA = new unsigned int[n];
	//unsigned int *SA;
	SA = (unsigned int*)malloc(sizeof(unsigned int)*n*2);
	if (s_ch == NULL || SA == NULL) {
		delete[] s_ch; delete[] SA;
		fprintf(stdout, "\nInsufficient memory, exit!");
		return;
	}


	// read the string into buffer.
	fprintf(stdout, "\nReading input string...");
	fseek(stdin, 0, SEEK_SET);
	for(i=0;i<n;i++) s_ch[i]=refer[i];
	// set the virtual sentinel
	s_ch[n - 1] = 0;


	auto start = std::chrono::high_resolution_clock::now();

    fprintf(stdout, "\nConstructing the suffix array...");
		auto SA_tmp = biovoltron::KISS1Sorter<uint32_t>{}.get_suffix_array(std::span{s_ch, n - 1});
		//SACA_K(s_ch, SA, n, 256, n, 0);
		fprintf(stdout, "\nsuccess kiss");
		// SACA_K(s_ch, SA, n, 256, n, 0);
		// fprintf(stdout, "\nsuccess SACA_K");
     

    printf("\nthe thread number is:%d", std::thread::hardware_concurrency());

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = finish - start;

    fprintf(stdout, "\nSize: %u bytes, Time: %5.3f seconds\n", n - 1, duration.count());
	/*
	 * SA_tmp is uint32_t type
	 * SA type is unsigend int
	 */

	for(int i = uint32_t{} ; i < n; i++)
		SA[i] = SA_tmp[i];
	SA[0] = n - 1;
	printf("sa_n=%d\n",n);
	*sa = SA;
}


unsigned int indenpendent_creadte_index(unsigned int text_length, char** input_refer, unsigned int compress_sa, char* filename)
{


	///����SAҪ��
	#define SA_counter_length 32

	///����Ҫ��
	///unsigned int compress_occ = 448, compress_SA_flag = 224;
	unsigned int compress_occ = 256, high_compress_occ=65536,compress_SA_flag = 224;
	typedef uint64_t bwt_string_type;

	///����SAҪ��
	///typedef uint32_t SA_flag_string_type;
	typedef uint64_t SA_flag_string_type;

	spdlog::info("compress_occ:{}", compress_occ);
	typedef uint32_t high_occ_table_type;
	unsigned int bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	unsigned int SA_flag_warp_number = sizeof(SA_flag_string_type)* 8;
	uint64_t* long_SA_flag = NULL;
	///����Ҫ��
	///unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned);
	unsigned int single_occ_bwt = sizeof(bwt_string_type) / sizeof(unsigned short);
	unsigned int occ_words = 4 / single_occ_bwt;
	unsigned int acctuall_bwt_gap = compress_occ + sizeof(unsigned int)* 4 * 8 / 2;

	///����SAҪ��
	///unsigned int acctuall_SA_flag_gap = compress_SA_flag + sizeof(SA_flag_string_type)* 8;
	unsigned int acctuall_SA_flag_gap = compress_SA_flag + SA_counter_length;
	///����SAҪ��
	unsigned int SA_counter_shift_length = sizeof(SA_flag_string_type)* 8 - SA_counter_length;

	spdlog::info("acctuall_SA_flag_gap:{}", acctuall_SA_flag_gap);
	bwt_string_type mode_4[4];
	bwt_string_type mode = (bwt_string_type)-1;
	bwt_string_type mode_high_1 = (bwt_string_type)1 << (SA_flag_warp_number - 1);
	bwt_string_type mode_16 = (bwt_string_type)65535;
	bwt_string_type mode_32 = ((bwt_string_type)-1) >> 32;
	bwt_string_type mode_high;
	bwt_string_type mode_low;
	SA_flag_string_type mode_SA_flag = (SA_flag_string_type)(((SA_flag_string_type)-1) << (sizeof(SA_flag_string_type)* 8 - 1));

	bwt_string_type pop_count_mode[4];
	unsigned int bwt_count_hash_table_bit = 16;
	unsigned int SA_length;
	unsigned int na, nc, ng, nt;
	unsigned int nacgt[4][258];
	unsigned int bwt_step = 1;
	unsigned int ctoi[256];

	unsigned int shapline;
	char filenames[200], filenameo[200], filenameb[200];
	char *refer = (*input_refer);

	bwt_string_type *bwt;
	unsigned int **occ;
	unsigned int *sa;

	char ch;
	unsigned int i, j;
	///����Ҫ��
	FILE *f1, *f2, *fs, *fb, *fo;
	ctoi['A'] = 0;
	ctoi['C'] = 1;
	ctoi['G'] = 2;
	ctoi['T'] = 3;
	ctoi['a'] = 0;
	ctoi['c'] = 1;
	ctoi['g'] = 2;
	ctoi['t'] = 3;



	unsigned int occ_line_number = 0;
	///this is the number of occ line
	occ_line_number = (text_length + 1) / compress_occ + 1;
	///this is the number of character which could be saved in one bwt_string_type
	bwt_warp_number = sizeof(bwt_string_type)* 8 / 2;
	///this is the number of bwt_string_type representing bwt string (does not include occ)
	unsigned int bwt_length = (text_length + 1) / bwt_warp_number + 1;
	///this is the number of bwt_string_type representing occ line
	///����Ҫ��
	/**
	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned int))
		/ (sizeof(bwt_string_type)) + 1;
	**/
	unsigned int occ_byte_length = (occ_line_number * 4 * sizeof(unsigned short))
		/ (sizeof(bwt_string_type)) + 1;



	bwt =
		(bwt_string_type *)malloc(sizeof(bwt_string_type)*(bwt_length + occ_byte_length + 1));





	///����SAҪ��
	/**
	unsigned int SA_flag_length = (text_length + 1) / SA_flag_warp_number + 1;
	unsigned int SA_flag_occ_length = (text_length + 1) / compress_SA_flag + 1;
	unsigned int SA_occ_byte_length = (SA_flag_occ_length * 1 * sizeof(unsigned int))
		/ (sizeof(SA_flag_string_type)) + 1;
	SA_flag_string_type* SA_flag =
		(SA_flag_string_type *)malloc(sizeof(SA_flag_string_type)*(SA_flag_length + SA_occ_byte_length + 1));
	**/
	unsigned int SA_flag_length = text_length + 1;
	unsigned int SA_flag_occ_length = (text_length + 1) / compress_SA_flag + 1;
	unsigned int SA_flag_byte_length = (SA_flag_length + SA_flag_occ_length * SA_counter_length + 1)
		/ SA_flag_warp_number + 1;
	SA_flag_string_type* SA_flag = (SA_flag_string_type *)malloc
		(sizeof(SA_flag_string_type)*SA_flag_byte_length);








	///����Ҫ��
	high_occ_table_type* high_occ_table;
	unsigned int high_occ_table_length = ((text_length + 1) / high_compress_occ + 1)*4+1;
	high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)*high_occ_table_length);

	///����SAҪ��
	///for (i = 0; i < SA_flag_length + SA_occ_byte_length + 1; i++)
	for (i = 0; i < SA_flag_byte_length; i++)
	{
		SA_flag[i] = (SA_flag_string_type)0;
	}





	for (i = 0; i < bwt_length + occ_byte_length + 1; i++)
	{
		bwt[i] = (bwt_string_type)0;
	}


	fprintf(stdout, "bwt_length=%llu,occ_byte_length=%llu, bwt_length+occ_byte_length=%llu\n",
		bwt_length, occ_byte_length, bwt_length + occ_byte_length);




	SA_length = text_length + 1;
	for (i = text_length; i > 0; i--)
	{
		ch = refer[i - 1];
		refer[i] = ctoi[ch] + 1;
	}
	refer[0] = 0;



	strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");
	f2 = fopen(filename, "w");
	fs = fopen(filenames, "w");
	fb = fopen(filenameb, "w");
	fo = fopen(filenameo, "w");




	indenpendent_get_sa_fromFILE(&sa, SA_length, &refer[1]);


	printf("SA has been generated!\n");
	na = nc = ng = nt = 0;
	for (i = 0; i<bwt_step; i++)
	for (j = 0; j <= (1 << (i + i + 2)); j++)
		nacgt[i][j] = 0;


	unsigned int bwt_iterater = 0;
	i = 0;
	unsigned int shift_length = 0;
	bwt_string_type tmp_bwt = (bwt_string_type)0;


	unsigned int tmp_occ = 0;
	///����Ҫ��
	/**
	unsigned int occ_bwt_number = (sizeof(unsigned int)* 4) / sizeof(bwt_string_type);
	**/
	unsigned int occ_bwt_number = (sizeof(unsigned short)* 4) / sizeof(bwt_string_type);
	for (bwt_iterater = 0; bwt_iterater < occ_bwt_number; bwt_iterater++)
	{
		bwt[bwt_iterater] = (bwt_string_type)0;
	}

	///����Ҫ��
	unsigned int high_occ_table_iterater = 0;
	for (high_occ_table_iterater = 0; high_occ_table_iterater < 4; high_occ_table_iterater++)
	{
		high_occ_table[high_occ_table_iterater] = 0;
	}




	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 258; j++)
		{
			nacgt[i][j] = 0;
		}
	}


	long long A_number[4] = { 0 };


	int flag = 0;

	i = 0;




	while (1)
	{

		tmp_bwt = (bwt_string_type)0;

		if (i >= SA_length)
		{
			break;
		}

		ch = refer[sa[i]] - 1;
		shift_length = (bwt_warp_number - i % bwt_warp_number - 1) * 2;
		tmp_bwt = (bwt_string_type)abs(ch);
		tmp_bwt = tmp_bwt << shift_length;





		bwt[bwt_iterater] = bwt[bwt_iterater] | tmp_bwt;





		tmp_bwt = (bwt_string_type)0;



		nacgt[0][ch + 1]++;



		A_number[abs(ch)]++;






		if (ch == -1 && i<SA_length)
		{
			shapline = i;
		}





		i++;



		if (i%bwt_warp_number == 0)
		{
			bwt_iterater++;
		}


		if (i % high_compress_occ == 0)
		{
			for (j = 0; j < 4; j++)
			{
				high_occ_table[high_occ_table_iterater] = nacgt[0][j + 1];
				high_occ_table_iterater++;
			}
		}


		if (i % compress_occ == 0)
		{


			///����Ҫ��,single_bwt_occ=4
			/**
			unsigned int single_bwt_occ = sizeof(bwt_string_type) / sizeof(unsigned int);
			**/
			unsigned int single_bwt_occ = sizeof(bwt_string_type) / sizeof(unsigned short);



			for (j = 0; j < 4; j++)
			{

				tmp_bwt = (bwt_string_type)0;

				///����Ҫ��
				///shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned int)* 8;
				///j=0, shift_length = 48
				///j=1, shift_length = 32
				///j=2, shift_length = 16
				///j=3, shift_length = 0
				shift_length = (single_bwt_occ - j % single_bwt_occ - 1) * sizeof(unsigned short)* 8;
				/**
				fprintf(stderr, "j=%llu\n", j);
				fprintf(stderr, "shift_length=%llu\n", shift_length);
				**/

				tmp_bwt = (bwt_string_type)nacgt[0][j + 1];

				/**
				fprintf(stderr, "tmp_bwt=%llu\n", tmp_bwt);
				**/

				///����Ҫ��
				tmp_bwt = tmp_bwt - high_occ_table[(i / high_compress_occ) * 4 + j];

				/**
				fprintf(stderr, "(i / high_compress_occ) * 4 + j=%llu\n", (i / high_compress_occ) * 4 + j);

				fprintf(stderr, "tmp_bwt=%llu\n", tmp_bwt);
				**/



				tmp_bwt = tmp_bwt << shift_length;





				bwt[bwt_iterater + j / single_bwt_occ] =
					bwt[bwt_iterater + j / single_bwt_occ] | tmp_bwt;


			}


			///����Ҫ��
			///bwt_iterater = bwt_iterater + single_bwt_occ;
			bwt_iterater = bwt_iterater + occ_words;
			/**
			fprintf(stderr, "single_bwt_occ=%llu\n", single_bwt_occ);
			fprintf(stderr, "occ_words=%llu\n", occ_words);
			**/




		}


		tmp_bwt = (bwt_string_type)0;


	}





	if (i%bwt_warp_number != 0 &&
		i % compress_occ != 0)
	{
		bwt_iterater++;
	}


	printf("Occ and BWT have been built!\n");



	fprintf(stdout, "write SA_length=%u, shapline=%u\n", SA_length, shapline);








	SA_flag_string_type tmp_SA_flag = (SA_flag_string_type)0;

	unsigned int SA_flag_iterater = 0;

	unsigned int SA_number = 0;

	///����SAҪ��
	unsigned int number_of_SA_flag_bits = 0;



	i = 0;

	SA_flag[SA_flag_iterater] = (SA_flag_string_type)0;


	///����SAҪ��
	///ע������SA_flag_iterater������++
	///SA_flag_iterater++;
	number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;

	/**
	while (1)
	{

		tmp_SA_flag = (SA_flag_string_type)0;

		if (i >= SA_length)
		{
			break;
		}



		if (sa[i] % compress_sa == 0)
		{
			tmp_SA_flag = (SA_flag_string_type)1;
			SA_number++;
		}
		else
		{
			tmp_SA_flag = (SA_flag_string_type)0;
		}


		shift_length = SA_flag_warp_number - i % SA_flag_warp_number - 1;

		tmp_SA_flag = tmp_SA_flag << shift_length;



		SA_flag[SA_flag_iterater] = SA_flag[SA_flag_iterater] | tmp_SA_flag;


		tmp_SA_flag = (SA_flag_string_type)0;



		i++;



		if (i%SA_flag_warp_number == 0)
		{
			SA_flag_iterater++;
		}


		if (i % compress_SA_flag == 0)
		{

			SA_flag[SA_flag_iterater] = SA_number;

			SA_flag_iterater++;

		}


		tmp_SA_flag = (SA_flag_string_type)0;


	}



	if (i%SA_flag_warp_number != 0 &&
		i % compress_SA_flag != 0)
	{
		SA_flag_iterater++;
	}
	**/


	while (1)
	{

		tmp_SA_flag = (SA_flag_string_type)0;

		if (i >= SA_length)
		{
			break;
		}



		if (sa[i] % compress_sa == 0)
		{
			tmp_SA_flag = (SA_flag_string_type)1;
			SA_number++;

		}
		else
		{
			tmp_SA_flag = (SA_flag_string_type)0;
		}



		///����SAҪ��
		///shift_length = SA_flag_warp_number - i % SA_flag_warp_number - 1;
		shift_length = SA_flag_warp_number - number_of_SA_flag_bits % SA_flag_warp_number - 1;

		tmp_SA_flag = tmp_SA_flag << shift_length;



		SA_flag[SA_flag_iterater] = SA_flag[SA_flag_iterater] | tmp_SA_flag;


		tmp_SA_flag = (SA_flag_string_type)0;



		///����SAҪ��
		number_of_SA_flag_bits++;

		i++;

		///����SAҪ��
		/**
		if (i%SA_flag_warp_number == 0)
		{
		SA_flag_iterater++;
		}


		if (i % compress_SA_flag == 0)
		{

		SA_flag[SA_flag_iterater] = SA_number;

		SA_flag_iterater++;

		}
		**/
		if (number_of_SA_flag_bits%SA_flag_warp_number == 0)
		{
			SA_flag_iterater++;
		}
		if (i % compress_SA_flag == 0)
		{

			if (number_of_SA_flag_bits % SA_flag_warp_number != 0)
			{
				fprintf(stderr, "ERROR! \n");
			}



			SA_flag[SA_flag_iterater] = SA_number;

			SA_flag[SA_flag_iterater] = (SA_flag_string_type)(SA_flag[SA_flag_iterater] << SA_counter_shift_length);

			number_of_SA_flag_bits = number_of_SA_flag_bits + SA_counter_length;

		}


		tmp_SA_flag = (SA_flag_string_type)0;


	}


	///����SAҪ��
	/**
	if (i%SA_flag_warp_number != 0 &&
	i % compress_SA_flag != 0)
	{
	SA_flag_iterater++;
	}
	**/
	if (number_of_SA_flag_bits%SA_flag_warp_number != 0)
	{
		SA_flag_iterater++;
	}

	///����SAҪ��,Ϊ�˷�ֹ�������ʱ���
	SA_flag_iterater++;



	///����Ҫ��
	fwrite(&high_occ_table_iterater, sizeof(unsigned int), 1, fo);
	fwrite(high_occ_table, sizeof(high_occ_table_type), high_occ_table_iterater, fo);
	printf("Occ has been writed!\n");



	fwrite(&bwt_iterater, sizeof(unsigned int), 1, fb);
	fwrite(bwt, sizeof(bwt_string_type), bwt_iterater, fb);
	printf("BWT has been writed!\n");

	fwrite(&SA_length, sizeof(SA_length), 1, f2);
	fwrite(&shapline, sizeof(shapline), 1, f2);

	printf("cc&sharp_line has been writed!\n");


	for (i = 0; i<bwt_step; i++)
	{

		nacgt[i][0] = 1;


		///fprintf(f2, "1");
		fwrite(&nacgt[i][0], sizeof(unsigned int), 1, f2);





		for (j = 1; j <= (1 << (i + i + 2)); j++)
		{
			nacgt[i][j] = nacgt[i][j] + nacgt[i][j - 1];

			fwrite(&nacgt[i][j], sizeof(unsigned int), 1, f2);
			fprintf(stdout, "nacgt[%u][%u]=%u\n", i, j, nacgt[i][j]);
			fflush(stdout);
		}
		///fprintf(f2, "\n");
	}


	fwrite(&SA_number, sizeof(SA_number), 1, fs);





	unsigned int ijkijkijkijk = 0;

	i = 0;

	unsigned int tmp_site = 0;






	while (1)
	{

		if (i >= SA_length)
		{
			break;
		}

		if (sa[i] % compress_sa == 0)
		{

			///if (compress_sa >= 4)
			///if (compress_sa >= 0)
			if (compress_sa >= 4)
			{
				///���ǰһλ��$����ô��sa��Ȼ��0����������ǲ�����Ҫ��
				//���Բ��ü�¼��ֻ��Ҫ�ж��ǲ���Ϊ0��Ϊ0�Ͷ���
				//�����ⶫ��Ӧ���Ǳ������010000000000000000000�����
				//tmp_site = (unsigned int)0;
				ch = refer[sa[i]] - 1;
				tmp_site = (unsigned int)abs(ch);

				tmp_site = tmp_site << 30;

				tmp_site = tmp_site | (sa[i] / compress_sa);


				fwrite(&tmp_site, sizeof(tmp_site), 1, fs);
			}
			else
			{
				fwrite(&sa[i], sizeof(sa[i]), 1, fs);
			}





			ijkijkijkijk++;
		}

		i++;

	}

	fprintf(stdout, "SA_number=%llu\n", SA_number);
	fprintf(stdout, "i=%llu, ijkijkijkijk=%llu\n", i, ijkijkijkijk);




	fprintf(stdout, "SA_flag_iterater=%llu\n", SA_flag_iterater);


	fwrite(&SA_flag_iterater, sizeof(SA_flag_iterater), 1, fs);
	fwrite(SA_flag, sizeof(SA_flag_string_type), SA_flag_iterater, fs);





	fwrite(&compress_sa, sizeof(unsigned int), 1, f2);
	fwrite(&compress_occ, sizeof(unsigned int), 1, f2);
	///����Ҫ��
	fwrite(&high_compress_occ, sizeof(unsigned int), 1, f2);







	fprintf(stdout, "Sucess!\n");
	fflush(stdout);



	fflush(fo);
	fflush(f2);
	fflush(fs);
	fflush(fb);

	fclose(fo);
	fclose(f2);
	fclose(fs);
	fclose(fb);
	printf("independent_creadte_index end!\n");
	return 0;


}



void init_bitmapper_index_params()
{
	unsigned int i;

	for (i = 0; i < 256; i++)
	{
		bitmapper_index_params.ctoi[i] = 4;
	}
	bitmapper_index_params.ctoi['A'] = 0;
	bitmapper_index_params.ctoi['C'] = 1;
	bitmapper_index_params.ctoi['G'] = 2;
	bitmapper_index_params.ctoi['T'] = 3;
	bitmapper_index_params.ctoi['a'] = 0;
	bitmapper_index_params.ctoi['c'] = 1;
	bitmapper_index_params.ctoi['g'] = 2;
	bitmapper_index_params.ctoi['t'] = 3;
	bitmapper_index_params.cut_thr = 1;

	bitmapper_index_params.tree_nodes_number = (pow(4, bitmapper_index_params.compress_sa) - 1) / 3;


	/**
	bitmapper_index_params.sp_tree = (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number);
	bitmapper_index_params.ep_tree = (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number);
	**/



	bitmapper_index_params.FMtree_queue
		= (unsigned int*)malloc(sizeof(unsigned int)* bitmapper_index_params.tree_nodes_number*3);



	bitmapper_index_params.FMtree_queue_start_point = 0;
	bitmapper_index_params.FMtree_queue_end_point = 0;




}







unsigned int load_index(char* filename_prefix)
{
	char filename[200], filename1[200], filenames[200], filenameo[200], filenameb[200];
	long long i, j, t;

	FILE *f1, *f2, *fs, *fb, *fout, *fo;



	bitmapper_index_params.pop_count_mode[0] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[0]
			= (bwt_string_type)(bitmapper_index_params.pop_count_mode[0] << 2) | (bwt_string_type)3;
	}


	bitmapper_index_params.pop_count_mode[1] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[1] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[1] << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.pop_count_mode[2] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[2] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[2] << 2) | (bwt_string_type)1;
	}

	bitmapper_index_params.pop_count_mode[3] = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.pop_count_mode[3] =
			(bwt_string_type)(bitmapper_index_params.pop_count_mode[3] << 2) | (bwt_string_type)0;
	}

	bitmapper_index_params.mode_high = (bwt_string_type)0;
	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_high =
			(bwt_string_type)(bitmapper_index_params.mode_high << 2) | (bwt_string_type)2;
	}


	bitmapper_index_params.mode_low = (bwt_string_type)0;

	for (i = 0; i < sizeof(bwt_string_type)* 8 / 2; i++)
	{
		bitmapper_index_params.mode_low =
			(bwt_string_type)(bitmapper_index_params.mode_low << 2) | (bwt_string_type)1;
	}



	char buffer[1024]; // 存储当前工作目录
	if (getcwd(buffer, sizeof(buffer)) != nullptr) {
		printf("当前工作目录:%s\n", buffer);
	} else {
		printf("获取当前工作目录失败:");
	}
	strcpy(filename, filename_prefix);
	//strcpy(filename + strlen(filename), ".index");
	strcpy(filenames, filename);
	strcpy(filenameo, filename);
	strcpy(filenameb, filename);
	strcpy(filenames + strlen(filename), ".sa");
	strcpy(filenameo + strlen(filename), ".occ");
	strcpy(filenameb + strlen(filename), ".bwt");


	///����Ҫ��
	fo = fopen(filenameo, "r");
	if (fo == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenameo);
		return 1;
	}

	fs = fopen(filenames, "r");
	if (fs == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenames);
		return 1;
	}
	fb = fopen(filenameb, "r");
	if (fb == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filenameb);
		return 1;
	}
	strcpy(filename1, filename);
	f2 = fopen(filename1, "r");

	if (f2 == NULL)
	{
		fprintf(stdout, "Failed to open %s! FMtree will exit ...\n", filename1);
		return 1;
	}


	fread(&bitmapper_index_params.SA_length, sizeof(bitmapper_index_params.SA_length), 1, f2);
	fread(&bitmapper_index_params.shapline, sizeof(bitmapper_index_params.shapline), 1, f2);

	fprintf(stdout, "shapline=%llu\n", bitmapper_index_params.shapline);



	for (j = 0; j <= 4; j++)
		fread(&bitmapper_index_params.nacgt[j], sizeof(bitmapper_index_params.nacgt[j]), 1, f2);




	fread(&bitmapper_index_params.compress_sa, sizeof(unsigned int), 1, f2);
	fread(&bitmapper_index_params.compress_occ, sizeof(unsigned int), 1, f2);
	fread(&bitmapper_index_params.high_compress_occ, sizeof(unsigned int), 1, f2);



	bitmapper_index_params.mode_4[0] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[1] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[2] = (bwt_string_type)-1;
	bitmapper_index_params.mode_4[3] = (bwt_string_type)0;


	fread(&bitmapper_index_params.bwt_length, sizeof(bitmapper_index_params.bwt_length), 1, fb);
	bitmapper_index_params.bwt = (bwt_string_type *)malloc(sizeof(bwt_string_type)*bitmapper_index_params.bwt_length);
	fread(bitmapper_index_params.bwt, sizeof(bwt_string_type), bitmapper_index_params.bwt_length, fb);
	printf("BWT has been loaded!\n");



	printf("SA_length=%u\n", bitmapper_index_params.SA_length);
	fread(&bitmapper_index_params.sparse_suffix_array_length,
		sizeof(bitmapper_index_params.sparse_suffix_array_length), 1, fs);
	printf("sparse_suffix_array_length=%u\n", bitmapper_index_params.sparse_suffix_array_length);
	bitmapper_index_params.sa
		= (unsigned int *)malloc(sizeof(unsigned int)*(bitmapper_index_params.sparse_suffix_array_length));
	fread(bitmapper_index_params.sa, sizeof(unsigned int), bitmapper_index_params.sparse_suffix_array_length, fs);








	fread(&bitmapper_index_params.SA_flag_iterater,
		sizeof(bitmapper_index_params.SA_flag_iterater), 1, fs);

	fprintf(stdout, "SA_flag_iterater=%llu\n", bitmapper_index_params.SA_flag_iterater);

	bitmapper_index_params.SA_flag =
		(SA_flag_string_type*)malloc(sizeof(SA_flag_string_type)*bitmapper_index_params.SA_flag_iterater);

	fread(bitmapper_index_params.SA_flag, sizeof(SA_flag_string_type),
		bitmapper_index_params.SA_flag_iterater, fs);


	fread(&bitmapper_index_params.high_occ_table_length, sizeof(bitmapper_index_params.high_occ_table_length), 1, fo);
	fprintf(stdout, "high_occ_table_length=%llu\n", bitmapper_index_params.high_occ_table_length);
	bitmapper_index_params.high_occ_table = (high_occ_table_type*)malloc(sizeof(high_occ_table_type)
		*bitmapper_index_params.high_occ_table_length);
	fread(bitmapper_index_params.high_occ_table, sizeof(high_occ_table_type),
		bitmapper_index_params.high_occ_table_length, fo);





	fclose(f2);
	fclose(fs);
	fclose(fb);
	fclose(fo);

	init_bitmapper_index_params();

	return 1;



}


























inline unsigned int bwt_get_sa_restrict_zero_steps_more_than_3
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;

	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
		*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;


	///����SAҪ��
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///����SAҪ��
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;




	///����SAҪ��
	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}


	///(*start) = sa[ans];
	(*start) = (sa[ans] & bitmapper_index_params.SA_header_mode)*bitmapper_index_params.compress_sa;


	return 1;

}




inline unsigned int bwt_get_sa_restrict_zero_steps
(SA_flag_string_type* SA_flag, unsigned int *sa, unsigned int sp, unsigned int* start)
{
	unsigned int l;
	unsigned int i = 0;
	///uint64_t tmp_SA_flag;
	unsigned int last;

	unsigned int actually_line;

	unsigned int ans;

	unsigned int last_wods;

	SA_flag_string_type tmp_SA_pop_count;


	l = sp;

	/**
	actually_line = ((l / bitmapper_index_params.compress_SA_flag)
		*bitmapper_index_params.acctuall_SA_flag_gap) / bitmapper_index_params.SA_flag_warp_number;
	**/


	actually_line = ((l / bitmapper_index_params.compress_SA_flag) << 8) >> 6;



	///����SAҪ��
	///last = l % compress_SA_flag;
	last = l % bitmapper_index_params.compress_SA_flag + SA_counter_length;

	///����SAҪ��
	///ans = SA_flag[actually_line];
	ans = SA_flag[actually_line] >> bitmapper_index_params.SA_counter_shift_length;


	///����SAҪ��
	if (last != 0)
	{
		///��Ҫlast�����Ѿ��ӹ�SA_counter_length�ˣ���������j=0�������Ӧ�ð�j = SA_counter_length
		unsigned int j = 0;

		tmp_SA_pop_count = SA_flag[actually_line] & bitmapper_index_params.SA_pop_count_mode;




		while (j + bitmapper_index_params.SA_flag_warp_number <= last)
		{

			ans = ans + __builtin_popcountll(tmp_SA_pop_count);

			j = j + bitmapper_index_params.SA_flag_warp_number;

			actually_line++;
			tmp_SA_pop_count = SA_flag[actually_line];

		}


		last = last - j;


		if (last != 0)
		{
			ans = ans +
				__builtin_popcountll((SA_flag_string_type)(tmp_SA_pop_count
				& (bitmapper_index_params.mode << (bitmapper_index_params.SA_flag_warp_number - last))));
		}


	}


	(*start) = sa[ans];

	return 1;

}


inline int bwt_accesss_pre_SA(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep,
	unsigned int* tmp_SA_length, unsigned int delta)
{
	unsigned int SA_start, SA_length, t, tmp_SA, tmp_SA_tail;

	if (ep - sp>1)
	{
		///����SAҪ��
		bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);

		/**
		if (SA_length>150000)
		{
			return 0;
		}
		**/


		for (t = 0; t<SA_length; t++)
		{
			tmp_SA = sa[t + SA_start];
			tmp_SA_tail = (tmp_SA & bitmapper_index_params.SA_header_mode)
				*bitmapper_index_params.compress_sa;



			if (tmp_SA_tail != 0 &&
				((bitmapper_index_params.SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)++] = tmp_SA_tail - 1;
			}


			///locates[(*tmp_SA_length)++] = sa[t + SA_start] & SA_header_mode - 1;

		}
	}
	else if (ep - sp == 1)   ///˵ʵ�������ǲ����ܵ������
	{




		///����SAҪ��
		if (bwt_get_sa_restrict_zero_steps(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{


			tmp_SA = locates[(*tmp_SA_length)];
			tmp_SA_tail = (tmp_SA & bitmapper_index_params.SA_header_mode)*bitmapper_index_params.compress_sa;



			if (tmp_SA_tail != 0 &&
				((bitmapper_index_params.SA_header_mode_reverse&tmp_SA) >> 30) == delta)
			{
				locates[(*tmp_SA_length)] = tmp_SA_tail - 1;
				(*tmp_SA_length)++;
			}



			///locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] & SA_header_mode - 1;
			///(*tmp_SA_length)++;

		}
	}

	return 1;

}



inline void bwt_accesss_SA_more_than_3(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;
	//printf("new locates\n");
	bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);
	// printf("SA_length=%u\n", SA_length);
	for (int t = 0; t<SA_length; t++)
	{
		locates[(*tmp_SA_length)++] = (sa[t + SA_start] & bitmapper_index_params.SA_header_mode)
			*bitmapper_index_params.compress_sa + occ;
	}
}


inline void bwt_accesss_SA_less_than_4(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	///����SAҪ��
	bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);

	for (t = 0; t<SA_length; t++)
	{
		locates[(*tmp_SA_length)++] = sa[t + SA_start] + occ;
	}


}


inline void bwt_accesss_SA_more_than_3_back_up(unsigned int *sa, SA_flag_string_type* SA_flag, unsigned int* locates,
	unsigned int sp, unsigned int ep, unsigned int occ,
	unsigned int* tmp_SA_length)
{
	unsigned int SA_start, SA_length, t;

	if (ep - sp>1)
	{


		///����SAҪ��
		bwt_direct_get_sa_interval_long(SA_flag, sp, ep, &SA_start, &SA_length);




		for (t = 0; t<SA_length; t++)
		{
			locates[(*tmp_SA_length)++] = (sa[t + SA_start] & bitmapper_index_params.SA_header_mode)
				*bitmapper_index_params.compress_sa + occ;
		}
	}
	else if (ep - sp == 1)
	{
		///����SAҪ��
		if (bwt_get_sa_restrict_zero_steps_more_than_3(SA_flag,
			sa, sp, (&(locates[(*tmp_SA_length)]))) == 1)
		{
			locates[(*tmp_SA_length)] = locates[(*tmp_SA_length)] + occ;
			(*tmp_SA_length)++;

		}

	}

}


inline void bwt_find_occ_all_sp_ep_optimal(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;
	unsigned int i_ep;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;

	unsigned int tmp_ans_C_sp, tmp_ans_T_sp;
	unsigned int tmp_ans_C_ep, tmp_ans_T_ep;
	unsigned int rank_1_sp, rank_1_ep;

	rank_1_sp = 0;
	rank_1_ep = 0;

	tmp_ans_C_sp = 0;
	tmp_ans_T_sp = 0;

	tmp_ans_C_ep = 0;
	tmp_ans_T_ep = 0;


	unsigned int high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];




	high_occ_table_line = (ep >> 16) << 2;
	(*ans_A_ep) = high_occ_table[high_occ_table_line];
	(*ans_C_ep) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_ep) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_ep) = high_occ_table[high_occ_table_line + 3];



	unsigned int actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;


	unsigned int actually_line_ep = ((ep >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;


	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}


	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;


	unsigned int need_line_sp = sp & 255;
	unsigned int need_line_ep = ep & 255;


	if (need_line_sp != 0)
	{

		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);

			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);



			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}


		///����ж�����ò����sp��ep������ͬһ��block
		///����ep��ֵò�ƿ���ֱ�Ӽ̳�
		if (flag == 1)
		{
			i_ep = i_sp;

			rank_1_ep = rank_1_sp;

			tmp_ans_C_ep = tmp_ans_C_sp;
		    (*ans_C_ep) = (*ans_C_sp);

			(*ans_G_ep) = (*ans_G_sp);

			tmp_ans_T_ep = tmp_ans_T_sp;
			(*ans_T_ep) = (*ans_T_sp);


			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;

			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);


			///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);


		}



		(*ans_C_sp) = (*ans_C_sp) + tmp_ans_C_sp;
		(*ans_G_sp) = (*ans_G_sp) + rank_1_sp - tmp_ans_C_sp - (tmp_ans_T_sp << 1);
		(*ans_T_sp) = (*ans_T_sp) + tmp_ans_T_sp;


		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];


















	if (need_line_ep != 0)
	{

		///����������ѭ������˵��sp��ep����һ��block
		///��ôֱ�ӽ���count�ͺ�, ���cache�����ʻ������ߵ�
		if (flag == 2)
		{



			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);

				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);


				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);

				P_tmp = P_tmp & P_B;


				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);




			}


			(*ans_C_ep) = (*ans_C_ep) + tmp_ans_C_ep;
			(*ans_G_ep) = (*ans_G_ep) + rank_1_ep - tmp_ans_C_ep - (tmp_ans_T_ep << 1);
			(*ans_T_ep) = (*ans_T_ep) + tmp_ans_T_ep;


			if ((bitmapper_index_params.shapline >=
				((ep >> 8) << 8))
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}












		}
		else
		{

			///�������sp��ep�ڲ�ͬ��block�������Ļ������һ��cache��ʧ
			(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

			(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

			(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

			(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

			actually_line_ep++;



			i_ep = 0;

			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);


				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);



				P_tmp = P_tmp & P_B;

				rank_1_ep = rank_1_ep + __builtin_popcountll(P_tmp);


				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_T_ep = tmp_ans_T_ep + __builtin_popcountll(P_A);



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				tmp_ans_C_ep = tmp_ans_C_ep + __builtin_popcountll(P_A);

			}


			(*ans_C_ep) = (*ans_C_ep) + tmp_ans_C_ep;
			(*ans_G_ep) = (*ans_G_ep) + rank_1_ep - tmp_ans_C_ep - (tmp_ans_T_ep<<1);
			(*ans_T_ep) = (*ans_T_ep) + tmp_ans_T_ep;



			if ((bitmapper_index_params.shapline >= ((ep >> 8) << 8)) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}

		}




	}
	else
	{
		(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

		(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

		(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

		(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

		actually_line_ep++;
	}



	(*ans_A_ep) += bitmapper_index_params.nacgt[0];
	(*ans_C_ep) += bitmapper_index_params.nacgt[1];
	(*ans_G_ep) += bitmapper_index_params.nacgt[2];
	(*ans_T_ep) += bitmapper_index_params.nacgt[3];








	return;
}





inline void bwt_find_occ_all_sp_ep_optimal_small(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;

	unsigned int tmp_ans_C_sp, tmp_ans_T_sp;
	unsigned int rank_1_sp = 0;


	tmp_ans_C_sp = 0;
	tmp_ans_T_sp = 0;


	unsigned int high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];


	unsigned int actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;



	///sp��block
	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;


	unsigned int need_line_sp = sp & 255;


	if (need_line_sp != 0)
	{

		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);

			///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);



			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}

		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;

			rank_1_sp = rank_1_sp + __builtin_popcountll(P_tmp);


			///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_T_sp = tmp_ans_T_sp + __builtin_popcountll(P_A);



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			tmp_ans_C_sp = tmp_ans_C_sp + __builtin_popcountll(P_A);


		}



		(*ans_C_sp) = (*ans_C_sp) + tmp_ans_C_sp;
		(*ans_G_sp) = (*ans_G_sp) + rank_1_sp - tmp_ans_C_sp - (tmp_ans_T_sp << 1);
		(*ans_T_sp) = (*ans_T_sp) + tmp_ans_T_sp;


		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];







	unsigned int delta;

	unsigned char count_ep[4];
	count_ep[0] = 0;
	count_ep[1] = 0;
	count_ep[2] = 0;
	count_ep[3] = 0;



	while (sp<ep)
	{
		actually_line_sp = ((sp >> 8)*bitmapper_index_params.acctuall_bwt_gap) + (sp & 255) + 32;


		delta = (bwt[(actually_line_sp >> 5)]
			>> ((bitmapper_index_params.bwt_warp_number - (actually_line_sp & 31) - 1) << 1))
			&(bwt_string_type)3;

		count_ep[delta]++;

		sp++;
	}


	(*ans_A_ep) = (*ans_A_sp) + count_ep[0];
	(*ans_C_ep) = (*ans_C_sp) + count_ep[1];
	(*ans_G_ep) = (*ans_G_sp) + count_ep[2];
	(*ans_T_ep) = (*ans_T_sp) + count_ep[3];

	return;
}




inline void bwt_find_occ_all_sp_ep_optimal_back_up(unsigned int sp, unsigned int ep, bwt_string_type *bwt, high_occ_table_type* high_occ_table,
	unsigned int* ans_A_sp, unsigned int* ans_C_sp, unsigned int* ans_G_sp, unsigned int* ans_T_sp,
	unsigned int* ans_A_ep, unsigned int* ans_C_ep, unsigned int* ans_G_ep, unsigned int* ans_T_ep)
{


	unsigned int i_sp;
	unsigned int i_ep;

	unsigned int flag = 0;

	bwt_string_type P_tmp, P_A, P_B;



	///unsigned int high_occ_table_line = (sp / bitmapper_index_params.high_compress_occ) * 4;
	///unsigned int high_occ_table_line = (sp / bitmapper_index_params.high_compress_occ) << 2;
	unsigned int high_occ_table_line = (sp >> 16) << 2;


	(*ans_A_sp) = high_occ_table[high_occ_table_line];
	(*ans_C_sp) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_sp) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_sp) = high_occ_table[high_occ_table_line + 3];



	///high_occ_table_line = (ep / bitmapper_index_params.high_compress_occ) * 4;
	///high_occ_table_line = (ep / bitmapper_index_params.high_compress_occ) << 2;
	high_occ_table_line = (ep >> 16) << 2;
	(*ans_A_ep) = high_occ_table[high_occ_table_line];
	(*ans_C_ep) = high_occ_table[high_occ_table_line + 1];
	(*ans_G_ep) = high_occ_table[high_occ_table_line + 2];
	(*ans_T_ep) = high_occ_table[high_occ_table_line + 3];
	///�Ӹ�����������ep��sp��superblock��һ��ĸ��ʻ��Ǻܴ��,����һ��ʼһ���cache�����ʻ��ǲ��͵�



	/**
	unsigned int actually_line_sp = ((sp / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) / bitmapper_index_params.bwt_warp_number;
	**/

	unsigned int actually_line_sp = ((sp >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;



	/**
	unsigned int actually_line_ep = ((ep / bitmapper_index_params.compress_occ)
		*bitmapper_index_params.acctuall_bwt_gap) / bitmapper_index_params.bwt_warp_number;
	**/

	unsigned int actually_line_ep = ((ep >> 8)
		*bitmapper_index_params.acctuall_bwt_gap) >> 5;










	if (actually_line_sp == actually_line_ep)
	{
		flag = 1;
	}






















	///sp��block
	(*ans_A_sp) += (bwt[actually_line_sp] >> 48);

	(*ans_C_sp) += (bwt[actually_line_sp] >> 32) & bitmapper_index_params.mode_16;

	(*ans_G_sp) += (bwt[actually_line_sp] >> 16) & bitmapper_index_params.mode_16;

	(*ans_T_sp) += bwt[actually_line_sp] & bitmapper_index_params.mode_16;

	actually_line_sp++;

	/**
	unsigned int need_line_sp = sp % bitmapper_index_params.compress_occ;
	unsigned int need_line_ep = ep % bitmapper_index_params.compress_occ;
	**/
	unsigned int need_line_sp = sp & 255;
	unsigned int need_line_ep = ep & 255;


	if (need_line_sp != 0)
	{




		i_sp = 0;

		while (i_sp + bitmapper_index_params.bwt_warp_number <= need_line_sp)
		{

			P_tmp = bwt[actually_line_sp++];


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);



			/**
			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);
			**/
			///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);





			i_sp = i_sp + bitmapper_index_params.bwt_warp_number;

		}


		///����ж�����ò����sp��ep������ͬһ��block
		///����ep��ֵò�ƿ���ֱ�Ӽ̳�
		if (flag == 1)
		{
			i_ep = i_sp;


			(*ans_C_ep) = (*ans_C_sp);

			(*ans_G_ep) = (*ans_G_sp);

			(*ans_T_ep) = (*ans_T_sp);

			actually_line_ep = actually_line_sp;

			flag = 2;
		}


		need_line_sp = need_line_sp - i_sp;

		if (need_line_sp != 0)
		{

			P_tmp = bwt[actually_line_sp++];
			/**
			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) * 2);
			**/

			P_B = bitmapper_index_params.mode
				<< ((bitmapper_index_params.bwt_warp_number - need_line_sp) << 1);


			P_tmp = P_tmp & P_B;



			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[1];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_C_sp) = (*ans_C_sp) + __builtin_popcountll(P_A);


			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[2];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_G_sp) = (*ans_G_sp) + __builtin_popcountll(P_A);

			/**
			P_A = P_tmp;
			P_B = P_A^bitmapper_index_params.pop_count_mode[3];
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);
			**/
			///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
			P_B = P_tmp;
			P_A = P_B >> 1;
			P_A = P_A & P_B;
			P_A = P_A & bitmapper_index_params.mode_low;
			(*ans_T_sp) = (*ans_T_sp) + __builtin_popcountll(P_A);






		}




		///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);



		/**
		if ((bitmapper_index_params.shapline >=
			(sp / bitmapper_index_params.compress_occ)*bitmapper_index_params.compress_occ)
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;
		**/


		if ((bitmapper_index_params.shapline >=
			((sp >> 8) << 8))
			&& (bitmapper_index_params.shapline<sp))
			(*ans_C_sp)--;


		(*ans_A_sp) = sp - (*ans_C_sp) - (*ans_G_sp) - (*ans_T_sp);


		if (sp>bitmapper_index_params.shapline)
		{
			(*ans_A_sp)--;
		}


	}



	(*ans_A_sp) += bitmapper_index_params.nacgt[0];
	(*ans_C_sp) += bitmapper_index_params.nacgt[1];
	(*ans_G_sp) += bitmapper_index_params.nacgt[2];
	(*ans_T_sp) += bitmapper_index_params.nacgt[3];


















	if (need_line_ep != 0)
	{

		///����������ѭ������˵��sp��ep����һ��block
		///��ôֱ�ӽ���count�ͺ�, ���cache�����ʻ������ߵ�
		if (flag == 2)
		{



			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);






				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];

				///P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) * 2);

				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);

				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);





			}




			///(*ans_A) = line - (*ans_C) - (*ans_G) - (*ans_T);



			/**
			if ((bitmapper_index_params.shapline >=
				(ep / bitmapper_index_params.compress_occ)*bitmapper_index_params.compress_occ)
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;
			**/
			if ((bitmapper_index_params.shapline >=
				((ep >> 8) << 8))
				&& (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}












		}
		else
		{

			///�������sp��ep�ڲ�ͬ��block�������Ļ������һ��cache��ʧ
			(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

			(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

			(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

			(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

			actually_line_ep++;



			i_ep = 0;

			while (i_ep + bitmapper_index_params.bwt_warp_number <= need_line_ep)
			{

				P_tmp = bwt[actually_line_ep++];


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);

				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);






				i_ep = i_ep + bitmapper_index_params.bwt_warp_number;

			}


			need_line_ep = need_line_ep - i_ep;

			if (need_line_ep != 0)
			{

				P_tmp = bwt[actually_line_ep++];


				///P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) * 2);
				P_B = bitmapper_index_params.mode << ((bitmapper_index_params.bwt_warp_number - need_line_ep) << 1);



				P_tmp = P_tmp & P_B;



				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[1];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_C_ep) = (*ans_C_ep) + __builtin_popcountll(P_A);


				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[2];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_G_ep) = (*ans_G_ep) + __builtin_popcountll(P_A);


				/**
				P_A = P_tmp;
				P_B = P_A^bitmapper_index_params.pop_count_mode[3];
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);
				**/
				///��Ϊ�κ������0����û�����κβ�������������������е㲻һ��������Ҫ���
				P_B = P_tmp;
				P_A = P_B >> 1;
				P_A = P_A & P_B;
				P_A = P_A & bitmapper_index_params.mode_low;
				(*ans_T_ep) = (*ans_T_ep) + __builtin_popcountll(P_A);





			}


			/**
			if ((bitmapper_index_params.shapline >= (ep / bitmapper_index_params.compress_occ)
				*bitmapper_index_params.compress_occ) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;
			**/
			if ((bitmapper_index_params.shapline >= ((ep >> 8)<<8)) && (bitmapper_index_params.shapline<ep))
				(*ans_C_ep)--;




			(*ans_A_ep) = ep - (*ans_C_ep) - (*ans_G_ep) - (*ans_T_ep);


			if (ep>bitmapper_index_params.shapline)
			{
				(*ans_A_ep)--;
			}

		}




	}
	else
	{
		(*ans_A_ep) += (bwt[actually_line_ep] >> 48);

		(*ans_C_ep) += (bwt[actually_line_ep] >> 32) & bitmapper_index_params.mode_16;

		(*ans_G_ep) += (bwt[actually_line_ep] >> 16) & bitmapper_index_params.mode_16;

		(*ans_T_ep) += bwt[actually_line_ep] & bitmapper_index_params.mode_16;

		actually_line_ep++;
	}



	(*ans_A_ep) += bitmapper_index_params.nacgt[0];
	(*ans_C_ep) += bitmapper_index_params.nacgt[1];
	(*ans_G_ep) += bitmapper_index_params.nacgt[2];
	(*ans_T_ep) += bitmapper_index_params.nacgt[3];








	return;
}


inline void enqueue_FMtree(unsigned int sp, unsigned int ep, unsigned int layer)
{
	unsigned int queue_point = bitmapper_index_params.FMtree_queue_end_point * 3;
	bitmapper_index_params.FMtree_queue[queue_point++] = sp;
	bitmapper_index_params.FMtree_queue[queue_point++] = ep;
	bitmapper_index_params.FMtree_queue[queue_point] = layer;


	bitmapper_index_params.FMtree_queue_end_point++;
}

inline void dequeue_FMtree(unsigned int* sp, unsigned int* ep, unsigned int* layer)

{

	unsigned int queue_point = bitmapper_index_params.FMtree_queue_start_point * 3;
	*sp =
		bitmapper_index_params.FMtree_queue[queue_point++];
	*ep =
		bitmapper_index_params.FMtree_queue[queue_point++];


	*layer =
		bitmapper_index_params.FMtree_queue[queue_point];


	bitmapper_index_params.FMtree_queue_start_point++;
}

inline unsigned int queue_length_FMtree()
{
	return (bitmapper_index_params.FMtree_queue_end_point - bitmapper_index_params.FMtree_queue_start_point);
}

inline void empty_queue_FMtree()
{
	bitmapper_index_params.FMtree_queue_start_point = 0;
	bitmapper_index_params.FMtree_queue_end_point = 0;
}


unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;

	unsigned int current_sp;
	unsigned int current_ep;
	unsigned int current_layer;

	unsigned int sp_A;
	unsigned int sp_C;
	unsigned int sp_G;
	unsigned int sp_T;


	unsigned int ep_A;
	unsigned int ep_C;
	unsigned int ep_G;
	unsigned int ep_T;



	///��Ϊ����early leaf node calculation
	unsigned int tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;

	///���ep-sp����ֵ֮�£�������λ��one-to-one���õ��������Ϳ���ֱ�ӷ�����
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		///ע�����������������չbitmapper_index_params.compress_sa - current_layer - 1
		///Ҳ����D-1��
		///�ͺ����ѭ����Ľ����һ��
		///�����ѭ������Ϊ�Ѿ�������early leaf node calculation, ����Ӧ����
		///tree_height - current_layer - 1
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			bitmapper_index_params.compress_sa - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///���������չ�����õ�sampled��λ��
	{
		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}


	///early leaf node calculation
	if (length_read >= 2)
	{
		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);



			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}


	}



	/**
	///early leaf node calculation
	if (length_read >= 2)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		///���ﲻҪ��
		///����SAҪ��
		if (bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta))
		{
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}
		}
		else
		{
			tree_height = bitmapper_index_params.compress_sa;
		}




	}
	**/





	///�����ÿ�
	empty_queue_FMtree();

	///������ִ�е����˵������Ҫ���ڵ�Ҫ��չ��
	///˵�����ĵ�һ������һ���Ѿ��������ˣ�����û���ҵ����е�λ��
	///���if�ж�˵���������������㼴D>2
	if (current_layer + 1 < tree_height)
	{
		/**
		if (current_ep - current_sp <= 4)
		{
			bwt_find_occ_all_sp_ep_optimal_small(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);
		}
		else
		{
			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);
		}
		**/


		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_A,
			&sp_C,
			&sp_G,
			&sp_T,
			&ep_A,
			&ep_C,
			&ep_G,
			&ep_T);





		if (ep_A != sp_A)
		{
			enqueue_FMtree(sp_A, ep_A, current_layer + 1);
		}

		if (ep_C != sp_C)
		{
			enqueue_FMtree(sp_C, ep_C, current_layer + 1);
		}

		if (ep_G != sp_G)
		{
			enqueue_FMtree(sp_G, ep_G, current_layer + 1);
		}
		if (ep_T != sp_T)
		{
			enqueue_FMtree(sp_T, ep_T, current_layer + 1);
		}
	}








	while (queue_length_FMtree()!=0)
	{

		////�ڵ����
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///���ep-sp����ֵ֮�£�������λ��one-to-one���õ�
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			///ע�����������������չtree_height - current_layer - 1
			///������bitmapper_index_params.compress_sa - current_layer - 1
			///Ҳ����D-2��
			///��ǰ���ѭ����Ľ����һ��
			///��Ϊ�����Ѿ�������early leaf node calculation, ����Ӧ����
			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///����������ֵ֮��,���䲻�����ٷ�֧��
			continue;




		}
		else  ///���������չ�����õ�sampled��λ��
		{

			bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{


			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);



			/**
			if (current_ep - current_sp <= 4)
			{
				bwt_find_occ_all_sp_ep_optimal_small(current_sp, current_ep,
					bitmapper_index_params.bwt,
					bitmapper_index_params.high_occ_table,
					&sp_A,
					&sp_C,
					&sp_G,
					&sp_T,
					&ep_A,
					&ep_C,
					&ep_G,
					&ep_T);
			}
			else
			{
				bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
					bitmapper_index_params.bwt,
					bitmapper_index_params.high_occ_table,
					&sp_A,
					&sp_C,
					&sp_G,
					&sp_T,
					&ep_A,
					&ep_C,
					&ep_G,
					&ep_T);
			}
			**/



			if (ep_A != sp_A)
			{

				enqueue_FMtree(sp_A, ep_A, current_layer + 1);
			}

			if (ep_C != sp_C)
			{
				enqueue_FMtree(sp_C, ep_C, current_layer + 1);
			}

			if (ep_G != sp_G)
			{
				enqueue_FMtree(sp_G, ep_G, current_layer + 1);
			}
			if (ep_T != sp_T)
			{
				enqueue_FMtree(sp_T, ep_T, current_layer + 1);
			}
		}



	}





	(*occurrences) = tmp_SA_length;
  return 1;
}


unsigned int locate_debug(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;

	unsigned int current_sp;
	unsigned int current_ep;
	unsigned int current_layer;

	unsigned int sp_A;
	unsigned int sp_C;
	unsigned int sp_G;
	unsigned int sp_T;


	unsigned int ep_A;
	unsigned int ep_C;
	unsigned int ep_G;
	unsigned int ep_T;



	///��Ϊ����early leaf node calculation
	unsigned int tree_height = bitmapper_index_params.compress_sa - 1;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;

	///���ep-sp����ֵ֮�£�������λ��one-to-one���õ��������Ϳ���ֱ�ӷ�����
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		///ע�����������������չbitmapper_index_params.compress_sa - current_layer - 1
		///Ҳ����D-1��
		///�ͺ����ѭ����Ľ����һ��
		///�����ѭ������Ϊ�Ѿ�������early leaf node calculation, ����Ӧ����
		///tree_height - current_layer - 1
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			bitmapper_index_params.compress_sa - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///���������չ�����õ�sampled��λ��
	{

		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}


	///early leaf node calculation
	if (length_read >= 2)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);



		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}


	}



	/**
	///early leaf node calculation
	if (length_read >= 2)
	{

	bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

	///���ﲻҪ��
	///����SAҪ��
	if (bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
	sp1, ep1,
	&tmp_SA_length, bitmapper_index_params.delta))
	{
	if (tmp_SA_length == number_of_hits)
	{
	(*occurrences) = tmp_SA_length;

	return 1;
	}
	}
	else
	{
	tree_height = bitmapper_index_params.compress_sa;
	}




	}
	**/





	///�����ÿ�
	empty_queue_FMtree();

	///������ִ�е����˵������Ҫ���ڵ�Ҫ��չ��
	///˵�����ĵ�һ������һ���Ѿ��������ˣ�����û���ҵ����е�λ��
	///���if�ж�˵���������������㼴D>2
	if (current_layer + 1 < tree_height)
	{
		debug_total++;
		if (current_ep - current_sp < 3
			&&
			current_ep - current_sp > 1)
		{
			debug_2++;
		}


		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_A,
			&sp_C,
			&sp_G,
			&sp_T,
			&ep_A,
			&ep_C,
			&ep_G,
			&ep_T);



		if (ep_A != sp_A)
		{
			enqueue_FMtree(sp_A, ep_A, current_layer + 1);
		}

		if (ep_C != sp_C)
		{
			enqueue_FMtree(sp_C, ep_C, current_layer + 1);
		}

		if (ep_G != sp_G)
		{
			enqueue_FMtree(sp_G, ep_G, current_layer + 1);
		}
		if (ep_T != sp_T)
		{
			enqueue_FMtree(sp_T, ep_T, current_layer + 1);
		}
	}








	while (queue_length_FMtree() != 0)
	{

		////�ڵ����
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///���ep-sp����ֵ֮�£�������λ��one-to-one���õ�
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			///ע�����������������չtree_height - current_layer - 1
			///������bitmapper_index_params.compress_sa - current_layer - 1
			///Ҳ����D-2��
			///��ǰ���ѭ����Ľ����һ��
			///��Ϊ�����Ѿ�������early leaf node calculation, ����Ӧ����
			bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///����������ֵ֮��,���䲻�����ٷ�֧��
			continue;




		}
		else  ///���������չ�����õ�sampled��λ��
		{

			bwt_accesss_SA_more_than_3(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{
			debug_total++;
			if (current_ep - current_sp < 3
				&&
				current_ep - current_sp > 1)
			{
				debug_2++;
			}

			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);


			if (ep_A != sp_A)
			{

				enqueue_FMtree(sp_A, ep_A, current_layer + 1);
			}

			if (ep_C != sp_C)
			{
				enqueue_FMtree(sp_C, ep_C, current_layer + 1);
			}

			if (ep_G != sp_G)
			{
				enqueue_FMtree(sp_G, ep_G, current_layer + 1);
			}
			if (ep_T != sp_T)
			{
				enqueue_FMtree(sp_T, ep_T, current_layer + 1);
			}
		}



	}





	(*occurrences) = tmp_SA_length;
  return 1;
}


unsigned int locate_less_than_4(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;

	unsigned int current_sp;
	unsigned int current_ep;
	unsigned int current_layer;

	unsigned int sp_A;
	unsigned int sp_C;
	unsigned int sp_G;
	unsigned int sp_T;


	unsigned int ep_A;
	unsigned int ep_C;
	unsigned int ep_G;
	unsigned int ep_T;



	///��Ϊ����early leaf node calculation
	unsigned int tree_height = bitmapper_index_params.compress_sa;


	current_sp = sp;
	current_ep = ep;
	current_layer = 0;

	///���ep-sp����ֵ֮�£�������λ��one-to-one���õ��������Ϳ���ֱ�ӷ�����
	if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
	{

		bwt_accesss_SA_cur_less_than_4(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length,
			tree_height - current_layer - 1);

		(*occurrences) = tmp_SA_length;

		return 1;


	}
	else  ///���������չ�����õ�sampled��λ��
	{

		bwt_accesss_SA_less_than_4(bitmapper_index_params.sa,
			bitmapper_index_params.SA_flag,
			locates,
			current_sp,
			current_ep,
			current_layer,
			&tmp_SA_length);

		///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
		if (tmp_SA_length == number_of_hits)
		{
			(*occurrences) = tmp_SA_length;

			return 1;
		}

	}





	///�����ÿ�
	empty_queue_FMtree();


	if (current_layer + 1 < tree_height)
	{
		bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
			bitmapper_index_params.bwt,
			bitmapper_index_params.high_occ_table,
			&sp_A,
			&sp_C,
			&sp_G,
			&sp_T,
			&ep_A,
			&ep_C,
			&ep_G,
			&ep_T);



		if (ep_A != sp_A)
		{
			enqueue_FMtree(sp_A, ep_A, current_layer + 1);
		}

		if (ep_C != sp_C)
		{
			enqueue_FMtree(sp_C, ep_C, current_layer + 1);
		}

		if (ep_G != sp_G)
		{
			enqueue_FMtree(sp_G, ep_G, current_layer + 1);
		}
		if (ep_T != sp_T)
		{
			enqueue_FMtree(sp_T, ep_T, current_layer + 1);
		}
	}








	while (queue_length_FMtree() != 0)
	{

		////�ڵ����
		dequeue_FMtree(&current_sp, &current_ep, &current_layer);

		///���ep-sp����ֵ֮�£�������λ��one-to-one���õ�
		if (current_ep - current_sp <= bitmapper_index_params.cut_thr)
		{

			bwt_accesss_SA_cur_less_than_4(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length,
				tree_height - current_layer - 1);


			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

			///����������ֵ֮��,���䲻�����ٷ�֧��
			continue;




		}
		else  ///���������չ�����õ�sampled��λ��
		{

			bwt_accesss_SA_less_than_4(bitmapper_index_params.sa,
				bitmapper_index_params.SA_flag,
				locates,
				current_sp,
				current_ep,
				current_layer,
				&tmp_SA_length);

			///����Ѿ��õ������е�λ�ã���ֱ�ӷ���
			if (tmp_SA_length == number_of_hits)
			{
				(*occurrences) = tmp_SA_length;

				return 1;
			}

		}




		if (current_layer + 1 < tree_height)
		{
			bwt_find_occ_all_sp_ep_optimal(current_sp, current_ep,
				bitmapper_index_params.bwt,
				bitmapper_index_params.high_occ_table,
				&sp_A,
				&sp_C,
				&sp_G,
				&sp_T,
				&ep_A,
				&ep_C,
				&ep_G,
				&ep_T);


			if (ep_A != sp_A)
			{

				enqueue_FMtree(sp_A, ep_A, current_layer + 1);
			}

			if (ep_C != sp_C)
			{
				enqueue_FMtree(sp_C, ep_C, current_layer + 1);
			}

			if (ep_G != sp_G)
			{
				enqueue_FMtree(sp_G, ep_G, current_layer + 1);
			}
			if (ep_T != sp_T)
			{
				enqueue_FMtree(sp_T, ep_T, current_layer + 1);
			}
		}



	}





	(*occurrences) = tmp_SA_length;
  return 1;
}


/**
unsigned int locate(char* pattern, unsigned int sp, unsigned int ep,
	unsigned int sp1, unsigned int ep1, unsigned int* locates, unsigned int length_read, unsigned int* occurrences)
{

	unsigned int tmp_SA_length = 0;
	unsigned int number_of_hits = ep - sp;
	unsigned int tree_layers_i = 0;
	unsigned int father_sp;
	unsigned int father_ep;

	unsigned int occ = 0;


	bitmapper_index_params.tree_index = 0;

	bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index] = sp;
	bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = ep;

	bitmapper_index_params.tree_layers = 0;

	bitmapper_index_params.need_step = bitmapper_index_params.compress_sa - bitmapper_index_params.tree_layers - 1;

	if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] -
		bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index] <= bitmapper_index_params.cut_thr)
	{
		///��������ڲ�Ҫ��
		///����SAҪ��
		bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
			locates,
			bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
			bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
			&tmp_SA_length, bitmapper_index_params.need_step);

		bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
			= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = (unsigned int)-1;
	}
	else
	{
		///���ﲻҪ��
		///����SAҪ��
		bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
			bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
			&tmp_SA_length);
	}




	if (length_read >= 2 &&
		tmp_SA_length != number_of_hits)
	{

		bitmapper_index_params.delta = bitmapper_index_params.ctoi[pattern[0]];

		///���ﲻҪ��
		///����SAҪ��
		bwt_accesss_pre_SA(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
			sp1, ep1,
			&tmp_SA_length, bitmapper_index_params.delta);


	}






	bitmapper_index_params.tree_layer_length = 1;







	bitmapper_index_params.tree_layer_length = 1;

	bitmapper_index_params.tree_index = 1;

	for (bitmapper_index_params.tree_layers = 1;
		bitmapper_index_params.tree_layers < bitmapper_index_params.compress_sa - 1;
		bitmapper_index_params.tree_layers++)
	{


		if (tmp_SA_length == number_of_hits)
		{
			break;
		}


		//need_step = compress_sa - tree_layers - 1;
		bitmapper_index_params.need_step
			= bitmapper_index_params.compress_sa - bitmapper_index_params.tree_layers - 2;

		occ = bitmapper_index_params.tree_layers;


		for (tree_layers_i = 0; tree_layers_i < bitmapper_index_params.tree_layer_length; tree_layers_i++)
		{

			father_sp = bitmapper_index_params.sp_tree
				[(bitmapper_index_params.tree_index - 1) / 4];
			father_ep = bitmapper_index_params.ep_tree
				[(bitmapper_index_params.tree_index - 1) / 4];





			if (father_sp == (unsigned int)-1)
			{
				bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
					= bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
					= (unsigned int)-1;
			}
			else
			{
				///����Ҫ��
				bwt_find_occ_all_sp_ep_optimal(father_sp, father_ep, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
					&bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2],
					&bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);





				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
					<= bitmapper_index_params.cut_thr)
				{
					///����Ҫ��
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);




					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]);

					fflush(stderr);



					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index] = (unsigned int)-1;
				}
				else
				{
					///���ﲻҪ��
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index], occ,
						&tmp_SA_length);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{

					break;
				}


				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
					<= bitmapper_index_params.cut_thr)
				{

					///����Ҫ��
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]);

					fflush(stderr);


					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1] = (unsigned int)-1;
				}
				else
				{
					///���ﲻҪ��
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1], occ,
						&tmp_SA_length);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 1],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 1]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{
					break;
				}


				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
					<= bitmapper_index_params.cut_thr)
				{

					///����Ҫ��
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);




					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]);

					fflush(stderr);




					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2] = (unsigned int)-1;
				}
				else
				{
					///���ﲻҪ��
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2], occ,
						&tmp_SA_length);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 2],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 2]);

					fflush(stderr);

				}




				if (tmp_SA_length == number_of_hits)
				{
					break;
				}

				if (bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]
					- bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
					<= bitmapper_index_params.cut_thr)
				{

					///����Ҫ��
					bwt_accesss_SA_cur_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, bitmapper_index_params.bwt, bitmapper_index_params.high_occ_table,
						locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3], occ,
						&tmp_SA_length, bitmapper_index_params.need_step);


					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);

					fflush(stderr);


					bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3]
						= bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]
						= (unsigned int)-1;
				}
				else
				{
					///���ﲻҪ��
					bwt_accesss_SA_more_than_3(bitmapper_index_params.sa, bitmapper_index_params.SA_flag, locates,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3], occ,
						&tmp_SA_length);



					fprintf(stderr, "tmp_SA_length=%llu, number_of_hits=%llu, current_layer=%llu, current_sp=%llu, current_ep=%llu\n",
						tmp_SA_length, number_of_hits, bitmapper_index_params.tree_layers,
						bitmapper_index_params.sp_tree[bitmapper_index_params.tree_index + 3],
						bitmapper_index_params.ep_tree[bitmapper_index_params.tree_index + 3]);

					fflush(stderr);

				}



				if (tmp_SA_length == number_of_hits)
				{
					break;
				}

			}

			bitmapper_index_params.tree_index = bitmapper_index_params.tree_index + 4;
		}

		bitmapper_index_params.tree_layer_length = bitmapper_index_params.tree_layer_length * 4;
	}

	(*occurrences) = tmp_SA_length;
}

**/

void debug_information()
{
	fprintf(stdout, "debug_2 = %llu \n", debug_2);
	fprintf(stdout, "debug_total = %llu \n", debug_total);
}
