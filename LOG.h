#pragma once
#define _DEFINE_DEPRECATED_HASH_CLASSES 0

#include <stack>
#include <vector>
#include <hash_map>
#include <boost/chrono.hpp>
#include <string>
#include <iostream>
#include <stdarg.h>
#include <Windows.h>
#include "MACROS.h"
#include "SCRIPT_READER.h"
#include "MULTITHREADING.h"
#include "FILE_UTILITIES.h"

using namespace stdext;

struct LOG_DATA
{
public: // Essential Data
	string									log_first_name;
	string									log_full_name;
	boost::chrono::system_clock::time_point begin_time_point;
	boost::chrono::system_clock::time_point end_time_point;
	boost::chrono::duration<double>			time;				// elapsed time for this frame

	boost::chrono::duration<double>			accumulated_time;	// accumulated time
	boost::chrono::duration<double>			average_time;

	boost::chrono::duration<double>			min_time;
	boost::chrono::duration<double>			max_time;
	int										num_time_accumulation;

public: // Constructor 
	LOG_DATA(const char* data_name)
		: log_full_name(data_name), num_time_accumulation(0), min_time((T)1000), max_time((T)0)
	{
		accumulated_time.zero();
	}
};

class LOG
{
public: // Essential Data
	// Members for evaluating accumulated & aveaged time
	typedef pair<string, LOG_DATA*>		log_pair;
	static hash_map<string, LOG_DATA*>	log_hash_map;

	// Ordinary members
	static stack<LOG_DATA*>				log_data_stack;
	
	static int							max_log_level;
	static int							current_log_level;

	static bool							use_log_file;
	static bool							use_stat;
	static bool							use_log_screen;
	static bool							display_accum_time;
	static bool							use_log;					// Whether file or screen log is turned on

	static char							va_buffer[1024];

	static string						log_path;
	static string						log_file_name;
	static string						stat_file_name;

	static string						level_tab;
	static ofstream						log_file_stream;
	static ofstream						stat_file_stream;

	static MULTITHREADING*				multithreading;

	static LOG							cout;

public: // Constructor and Destructor
	LOG(void)
	{}

	~LOG(void)
	{
		if (log_hash_map.size() != 0)
		{
			hash_map<string, LOG_DATA*>::iterator iter_hash_map = log_hash_map.begin();
			for (int i = 0; i < log_data_stack.size(); i++)
			{
				LOG_DATA* log_data = iter_hash_map->second;
				DELETE_POINTER(log_data);
				iter_hash_map++;
			}
			log_hash_map.clear();
		}

		if (log_data_stack.size() != 0)
		{
			for (int i = 0; i < log_data_stack.size(); i++)
			{
				log_data_stack.pop();
			}
		}

		log_file_stream.close();
		stat_file_stream.close();
	}

public: // Initialization Functions
	static void Initialize(MULTITHREADING* multithreading_input, const bool& out_screen, const bool& out_file, const char* file_path, const int& max_level)
	{
		multithreading = multithreading_input;

		max_log_level = max_level;

		use_log_file = out_file;
		use_log_screen = out_screen;

		log_file_name = file_path;

		if (use_log_file == true)
		{
			log_file_stream.open(log_file_name.data(), ios_base::trunc || ios_base::out);
			assert(log_file_stream.is_open() == true);
		}

		if (use_log_file == true || use_log_screen == true)
		{
			use_log = true;
		}

		else
		{
			use_log = false;
		}
	}

	static void InitializeFromScriptBlock(MULTITHREADING* multithreading_input, const SCRIPT_BLOCK& outer_block, const string script_file_dir)
	{
		multithreading = multithreading_input;

		max_log_level		= outer_block.GetInteger("max_log_level", 10000);
		use_log_file		= outer_block.GetBoolean("use_log_file", false);
		use_stat			= outer_block.GetBoolean("use_stat", false);
		use_log_screen		= outer_block.GetBoolean("use_log_screen", false);
		display_accum_time	= outer_block.GetBoolean("display_accum_time", false);

		// When the common path is not set
		log_path = script_file_dir;
		log_path.append(outer_block.GetString("output_log_path", "logs"));
		log_file_name = outer_block.GetString("output_log_file_name", "log");

		FILE_UTILITIES file_utilites;
		if (!file_utilites.Directory_Exists(log_path))
		{
			file_utilites.Create_Directory(log_path);
		}

		if (use_log_file == true)
		{
			string log_file_dest;
			log_file_dest.append(log_path);
			log_file_dest.append("\\");
			log_file_dest.append(log_file_name);
			log_file_dest.append(".log");

			log_file_stream.open(log_file_dest.data(), ios_base::trunc || ios_base::out);
			assert(log_file_stream.is_open() == true);
		}

		if (use_stat == true)
		{
			string stat_file_dest;
			stat_file_dest.append(log_path);
			stat_file_dest.append("\\");
			stat_file_dest.append(log_file_name);
			stat_file_dest.append(".stat");

			stat_file_stream.open(stat_file_dest.data(), ios_base::trunc || ios_base::out);
			assert(stat_file_stream.is_open() == true);
		}

		if (use_log_file == true || use_log_screen == true)
		{
			use_log = true;
		}
		
		else
		{
			use_log = false;
		}
	}

	static void Begin(const char* log_format, ...)
	{
		static int cnt = 0;
		if (use_log == false)
		{
			return;
		}
		
		if (++current_log_level > max_log_level)
		{
			return;
		}

		va_list ap;
		va_start(ap, log_format);
			vsprintf(va_buffer, log_format, ap);
		va_end(ap);

		// 1. Saving log name
		string log_full_name = va_buffer;
		char* first_word = strtok(va_buffer, " ");
		// Tp save following ones such as "frame 1", "frame 2" be same label
		string log_first_name(first_word);

		// 2. Using both log_hash_map & log_data_stack
		if (display_accum_time)
		{
			hash_map<string, LOG_DATA*>::iterator iter_hash_map;
			iter_hash_map = log_hash_map.find(log_first_name);

			if (iter_hash_map == log_hash_map.end())		// 2a. new log : inserting
			{
				// LOG_DATA making
				LOG_DATA* log_data = new LOG_DATA(va_buffer);
				log_data->begin_time_point = boost::chrono::system_clock::now();
				log_data->log_first_name = log_first_name;
				log_data->log_full_name = log_full_name;

				log_hash_map.insert(log_pair(log_first_name, log_data));
				log_data_stack.push(log_data);
			}
			else                                            // 2b. existing log : retrieve
			{
				LOG_DATA* log_data = iter_hash_map->second;
				log_data->begin_time_point = boost::chrono::system_clock::now();
				log_data->log_full_name = log_full_name;

				log_data_stack.push(iter_hash_map->second);
			}
		}

		// 3. Only using log_data_stack
		else
		{
			// LOG DATA making
			LOG_DATA* log_data = new LOG_DATA(va_buffer);
			log_data->begin_time_point = boost::chrono::system_clock::now();
			log_data->log_first_name = log_first_name;
			log_data->log_full_name = log_full_name;

			// 4. LOG_DATA stacking
			log_data_stack.push(log_data);
		}

		// 4a. Screen cout
		if (use_log_screen == true)
		{
			std::cout << level_tab << "|BEGIN| " << log_full_name << std::endl;
		}
		
		// 4b. File out
		if (use_log_file == true)
		{
			log_file_stream << level_tab << "|BEGIN| " << log_full_name << std::endl;
		}

		level_tab.append("  ");
	}
		
	static void End()
	{
		if (use_log == false)
		{
			return;
		}
			
		if (current_log_level-- > max_log_level)
		{
			return;
		}
		assert(current_log_level >= -1);

		level_tab.erase(0, 2);
		
		LOG_DATA* log_data = log_data_stack.top();
		log_data_stack.pop();
		log_data->end_time_point = boost::chrono::system_clock::now();
		log_data->time           = log_data->end_time_point - log_data->begin_time_point;

		if (log_data->min_time > log_data->time)
		{
			log_data->min_time = log_data->time;
		}
		if (log_data->max_time < log_data->time)
		{
			log_data->max_time = log_data->time;
		}

		if (display_accum_time)
		{
			// 1. log time information updata
			hash_map<string, LOG_DATA*>::iterator iter_hash_map = log_hash_map.find(log_data->log_first_name);
			if (iter_hash_map == log_hash_map.end())
			{
				std::cout << "LOG error!!!" << std::endl;
			}

			LOG_DATA* log_data = iter_hash_map->second;
			log_data->accumulated_time += log_data->time;
			log_data->num_time_accumulation++;
			log_data->average_time = log_data->accumulated_time/(double)log_data->num_time_accumulation;

			// 2. output log
			if (use_log_screen == true)
			{
				std::cout << level_tab << "|END|	" << log_data->log_full_name << " [this:" << log_data->time.count() << "] [ave: " << log_data->average_time.count() << "] [accum:" << log_data->accumulated_time.count() << "]";
				std::cout << "[min:" << log_data->min_time.count() << "] [max:" << log_data->max_time.count() << "]" << std::endl;
			}
			if (use_log_file == true)
			{
				log_file_stream << level_tab << "|END|  " << log_data->log_full_name << " [this:" << log_data->time.count() << "] [ave:" << log_data->average_time.count() << "] [accum:" << log_data->accumulated_time.count() << "]";
				log_file_stream << "[min:" << log_data->min_time.count() << "] [max:" << log_data->max_time.count() << "]" << std::endl;
			}
		}
		else
		{
			// 2. Output log
			if (use_log_screen == true)
			{
				std::cout << level_tab << "|END|  " << log_data->log_full_name << " [" << log_data->time.count() << "]" << std::endl;
			}
			if (use_log_file == true)
			{
				log_file_stream << level_tab << "|END|  " << log_data->log_full_name << " [" << log_data->time.count() << "]" << std::endl;
			}
			DELETE_POINTER(log_data);
		}
	}

	static void Begin(int thread_id, const char* log_format, ...)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			va_list ap;
			va_start(ap, log_format);
				vsprintf(va_buffer, log_format, ap);
			va_end(ap);

			Begin(va_buffer);
		}
		END_HEAD_THREAD_WORK
	}

	static void End(int thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			End();
		}
		END_HEAD_THREAD_WORK
	}

	static void Stat(const char* log_format, ...)
	{
		if (use_stat == false)
		{
			return;
		}

		va_list ap;
		va_start(ap, log_format);
			vsprintf(va_buffer, log_format, ap);
		va_end(ap);

		stat_file_stream << va_buffer << std::endl;
	}

	static void Stat(int thread_id, const char* log_format, ...)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			va_list ap;
			va_start(ap, log_format);
				vsprintf(va_buffer, log_format, ap);
			va_end(ap);

			Stat(va_buffer);
		}
		END_HEAD_THREAD_WORK
	}

	static void Begin(boost::chrono::system_clock::time_point& start_time_input, const char* log_name_input)
	{
		start_time_input = boost::chrono::system_clock::now();

		if (use_log_screen == true)
		{
			std::cout << "|BEGIN| " << log_name_input << std::endl;
		}
		if (use_log_file == true)
		{
			log_file_stream << "|BEGIN| " << log_name_input << std::endl;
		}
	}

	static void End(boost::chrono::system_clock::time_point start_time_input, T& elapsed_time_input, const char* log_name_input, T& accum_time_input, int frame_input)
	{
		boost::chrono::duration<double> elapsed_time;
		elapsed_time = boost::chrono::system_clock::now() - start_time_input;
		elapsed_time_input = (T)elapsed_time.count();
		accum_time_input += elapsed_time_input;
		T ave_time_input = accum_time_input/(T)frame_input;

		if (use_log_screen == true)
		{
			std::cout << "|END|  " << log_name_input << " [this:" << elapsed_time_input << "] [ave:" << ave_time_input << "]" << "] [accum:" << accum_time_input << "] frame:" << frame_input << std::endl;
		}
		if (use_log_file == true)
		{
			log_file_stream << "|END|  " << log_name_input << " [this:" << elapsed_time_input << "] [ave:" << ave_time_input << "] [accum:" << accum_time_input << "] frame:" << frame_input << std::endl;
		}
	}

	static void Begin(int thread_id, boost::chrono::system_clock::time_point& start_time_input, const char* log_name_input)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			start_time_input = boost::chrono::system_clock::now();

			if (use_log_screen == true)
			{
				std::cout << "|BEGIN| " << log_name_input << std::endl;
			}
			if (use_log_file == true)
			{
				log_file_stream << "|BEGIN| " << log_name_input << std::endl;
			}
		}
		END_HEAD_THREAD_WORK
	}

	static void End(int thread_id, boost::chrono::system_clock::time_point start_time_input, T& elapsed_time_input, const char* log_name_input, T& accum_time_input, int frame_input)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			boost::chrono::duration<double> elapsed_time;
			elapsed_time = boost::chrono::system_clock::now() - start_time_input;
			elapsed_time_input = (T)elapsed_time.count();
			accum_time_input += elapsed_time_input;
			T ave_time_input = accum_time_input/(T)frame_input;

			if (use_log_screen == true)
			{
				std::cout << "|END|  " << log_name_input << " [this:" << elapsed_time_input << "][ave:" << ave_time_input << "]" << "][accum:" << accum_time_input << "] frame:" << frame_input << std::endl;
			}
			if (use_log_file == true)
			{
				log_file_stream << "|END|  " << log_name_input << " [this:" << elapsed_time_input << "][ave:" << ave_time_input << "][accum:" << accum_time_input << "] frame:" << frame_input << std::endl;
			}
		}
		END_HEAD_THREAD_WORK
	}
};

LOG& operator<< (LOG& out, char c);
LOG& operator<< (LOG& out, signed char c);
LOG& operator<< (LOG& out, unsigned char c);

LOG& operator<< (LOG& out, const char* s);
LOG& operator<< (LOG& out, const signed char* s);
LOG& operator<< (LOG& out, const unsigned char* s);

LOG& operator<< (LOG& out, bool val);
LOG& operator<< (LOG& out, short val);
LOG& operator<< (LOG& out, unsigned short val);
LOG& operator<< (LOG& out, int val);
LOG& operator<< (LOG& out, unsigned int val);
LOG& operator<< (LOG& out, long val);
LOG& operator<< (LOG& out, unsigned long val);
LOG& operator<< (LOG& out, float val);
LOG& operator<< (LOG& out, double val);
LOG& operator<< (LOG& out, long double val);
LOG& operator<< (LOG& out, const void* val);
LOG& operator<< (LOG& out, ostream& (*pf)(ostream&));
LOG& operator<< (LOG& out, ios& (*pf)(ios&));
LOG& operator<< (LOG& out, ios_base& (*pf)(ios_base&));