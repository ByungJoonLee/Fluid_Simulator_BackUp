#include "stdafx.h"
#include "LOG.h"

hash_map<string, LOG_DATA*> LOG::log_hash_map;
stack<LOG_DATA*> LOG::log_data_stack;

int LOG::max_log_level = 10000;
int LOG::current_log_level = -1;

bool LOG::use_log_file = false;
bool LOG::use_log_screen = false;
bool LOG::display_accum_time = false;
bool LOG::use_log = false;
bool LOG::use_stat = false;

char LOG::va_buffer[1024];

string LOG::log_path;
string LOG::log_file_name;

string LOG::level_tab;

ofstream LOG::log_file_stream;
ofstream LOG::stat_file_stream;

MULTITHREADING* LOG::multithreading = 0;

LOG LOG::cout;

LOG& operator<< (LOG& out, char c)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << c;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << c;
	}
	return out;
}

LOG& operator<< (LOG& out, signed char c)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << c;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << c;
	}
	return out;
}

LOG& operator<< (LOG& out, unsigned char c)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << c;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << c;
	}
	return out;
}

LOG& operator<< (LOG& out, const char* s)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << s;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << s;
	}
	return out;
}

LOG& operator<< (LOG& out, const signed char* s)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << s;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << s;
	}
	return out;
}

LOG& operator<< (LOG& out, const unsigned char* s)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << s;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << s;
	}
	return out;
}

LOG& operator<< (LOG& out, bool val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, short val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, unsigned short val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}	

LOG& operator<< (LOG& out, int val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, unsigned int val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, long val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, unsigned long val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, float val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, double val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, long double val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, const void* val)
{
	if (LOG::use_log_screen == true)
	{
		std::cout << val;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << val;
	}
	return out;
}

LOG& operator<< (LOG& out, ostream& (*pf)(ostream&))
{
	if (LOG::use_log_screen == true)
	{
		std::cout << pf;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << pf;
	}
	return out;
}

LOG& operator<< (LOG& out, ios& (*pf)(ios&))
{
	if (LOG::use_log_screen == true)
	{
		std::cout << pf;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << pf;
	}
	return out;
}

LOG& operator<< (LOG& out, ios_base& (*pf)(ios_base&))
{
	if (LOG::use_log_screen == true)
	{
		std::cout << pf;
	}
	if (LOG::use_log_file == true)
	{
		LOG::log_file_stream << pf;
	}
	return out;
}