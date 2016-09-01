#include "function_store.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include "npn_canonization.hpp"
#include <boost/format.hpp>


using namespace std;
using namespace cirkit;

namespace majesty {

	string function_store::compute_min_size_depth_xmg(const tt& function) {
		return compute_min_size_depth_xmg(function, ALL);
	}

	string function_store::compute_min_size_depth_xmg(const tt& function, NODE_TYPE type) {
		auto command = "cirkit -ec \"tt -l " + to_string(function) + ";"; 
		if (type == MAJ) {
			command += "exact_mig -v; ps -x; convert --mig_to_expr; print -e\"";
		} else {
			command += "exact_xmg -v; ps -x; convert --xmg_to_expr; print -e\"";
		}
#ifndef _WIN32
		command += " > " + _logfile;
#endif
		auto retval = system(command.c_str());
		if (retval != 0) {
			throw "Exact synthesis through Cirkit failed";
		}
		_remove_log = true;

		string line, lastline;
		ifstream filein;
		filein.open(_logfile);
		while (getline(filein, line)) {
			lastline = line;
		}
		filein.close();
		boost::algorithm::trim_right(lastline);

		return lastline;
	}
	
	function_store::~function_store() {
#ifndef _WIN32
		if (_rcontext != NULL) {
			redisFree(_rcontext);
		}
		if (_remove_log) {
		   	if (remove(_logfile.c_str())) {
				cerr << "Warning: unable to remove temporary file!" << endl;
			}
		}
#endif
	}

	function_store::function_store() : function_store("127.0.0.1", 6379) {
	}

	function_store::function_store(const std::string& server_url, 
			const unsigned port) {
		auto path = boost::filesystem::unique_path();
		_logfile = path.native();
		cout << "connecting to " << server_url << ":" << port << endl;
#ifndef _WIN32
		_rcontext = redisConnect(server_url.c_str(), port);
		if (_rcontext == NULL || _rcontext->err) {
			if (_rcontext) {
				throw "Error initializing Redis context";
			} else {
				throw "Error allocating Redis context";
			}
			exit(1);
		}
#endif
	}

	std::string function_store::min_size_depth_xmg(const cirkit::tt& f) {
		return min_size_depth_xmg(f, ALL);
	}

	std::string function_store::min_size_depth_xmg(const cirkit::tt& f, NODE_TYPE type) {
		// First see if the optimum xmg has already been computed
		string expr;
#ifndef _WIN32
		redisReply* reply = (redisReply*)
			redisCommand(_rcontext, "GET %s:expr", to_string(f).c_str());
		assert(reply != NULL); // NOTE: handle errors!
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				expr = compute_min_size_depth_xmg(f, type);
				reply = (redisReply*)redisCommand(_rcontext, "SET %s:expr %s", 
						to_string(f).c_str(), expr.c_str());
				assert(reply != NULL); // NOTE: handle errors!
				break;
			case REDIS_REPLY_STRING:
				expr = string(reply->str, reply->len);
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				throw errorstring;
				break;
		}
		freeReplyObject(reply);
#endif
		return expr;
	}
	
	tt function_store::npn_canon(const tt& f, tt& phase, 
			vector<unsigned>& perm) {
		tt npn;
#ifndef _WIN32
		redisReply* reply = (redisReply*)
			redisCommand(_rcontext, "GET %s:npn", to_string(f).c_str());
		assert(reply != NULL); // NOTE: handle errors!
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				npn = exact_npn_canonization(f, phase, perm);
				assert(perm.size() == tt_num_vars(f));
				reply = (redisReply*)redisCommand(_rcontext, "MULTI");
				assert(reply != NULL && reply->type == REDIS_REPLY_STATUS);
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"SET %s:npn %s", to_string(f).c_str(), 
						to_string(npn).c_str());
				assert(reply != NULL);
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"SET %s:phase %s", to_string(f).c_str(), 
						to_string(phase).c_str());
				assert(reply != NULL);
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, "DEL %s:perm",
						to_string(f).c_str());
				assert(reply != NULL);
				freeReplyObject(reply);
				for (auto i = 0u; i < perm.size(); i++) {
					reply = (redisReply*)redisCommand(_rcontext, 
							"RPUSH %s:perm %u", to_string(f).c_str(), perm[i]);
					assert(reply != NULL);
					freeReplyObject(reply);
				}
				reply = (redisReply*)redisCommand(_rcontext, "EXEC");
				assert(reply != NULL);
				break;
			case REDIS_REPLY_STRING:
				npn = tt(string(reply->str, reply->len));
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"GET %s:phase", to_string(f).c_str());
				assert(reply != NULL);
				phase = tt(string(reply->str, reply->len));
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"LRANGE %s:perm 0 -1", to_string(f).c_str());
				assert(reply != NULL && reply->type == REDIS_REPLY_ARRAY);
				for (auto i = 0u; i < reply->elements; i++) {
					auto preply = reply->element[i];
					auto permstring = string(preply->str, preply->len);
					perm.push_back(stoi(permstring));
				}
				assert(perm.size() == tt_num_vars(f));
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				throw errorstring;
				break;
		}
		freeReplyObject(reply);
#endif
		return npn;
	}

}
