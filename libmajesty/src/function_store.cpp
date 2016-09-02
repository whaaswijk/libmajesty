#include "function_store.h"
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <chrono>
#include "npn_canonization.hpp"
#include <boost/format.hpp>
#include <convert.h>


using namespace std;
using namespace cirkit;
using boost::optional;

namespace majesty {

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
		cout << "context: " << _rcontext << endl;
		cout << "context->err: " << _rcontext->err << endl;
		if (_rcontext == NULL || _rcontext->err) {
			if (_rcontext) {
				throw runtime_error("Error initializing Redis context");
			} else {
				throw runtime_error("Error allocating Redis context");
			}
		}
#endif
	}

	string function_store::min_size_xmg(const cirkit::tt& f) {
		return min_size_xmg(f, 0).value();
	}

	optional<string> function_store::min_size_xmg(const cirkit::tt& f, unsigned timeout) {
		// First see if the optimum xmg has already been computed
		optional<string> expr;
#ifndef _WIN32
		redisReply* reply = (redisReply*)
			redisCommand(_rcontext, "GET %s:expr", to_string(f).c_str());
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				expr = exact_xmg_expression(f, timeout);
				if (expr) { // A timeout may have occurred on synthesis
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:expr %s",
						to_string(f).c_str(), expr.c_str());
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
				}
				break;
			case REDIS_REPLY_STRING:
				expr = string(reply->str, reply->len);
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				throw runtime_error(errorstring);
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
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				npn = exact_npn_canonization(f, phase, perm);
                if (perm.size() != tt_num_vars(f)) {
                    throw runtime_error("Error: perm.size() != tt_num_vars(f)");
                }
				reply = (redisReply*)redisCommand(_rcontext, "MULTI");
                if (reply == NULL || reply->type != REDIS_REPLY_STATUS) {
                    throw runtime_error("Error connecting to server");
                }
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"SET %s:npn %s", to_string(f).c_str(), 
						to_string(npn).c_str());
                if (reply == NULL) {
                    throw runtime_error("Error connecting to server");
                }
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"SET %s:phase %s", to_string(f).c_str(), 
						to_string(phase).c_str());
                if (reply == NULL) {
                    throw runtime_error("Error connecting to server");
                }
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, "DEL %s:perm",
						to_string(f).c_str());
                if (reply == NULL) {
                    throw runtime_error("Error connecting to server");
                }
				freeReplyObject(reply);
				for (auto i = 0u; i < perm.size(); i++) {
					reply = (redisReply*)redisCommand(_rcontext, 
							"RPUSH %s:perm %u", to_string(f).c_str(), perm[i]);
                    if (reply == NULL) {
                        throw runtime_error("Error connecting to server");
                    }
					freeReplyObject(reply);
				}
				reply = (redisReply*)redisCommand(_rcontext, "EXEC");
                if (reply == NULL) {
                    throw runtime_error("Error connecting to server");
                }
				break;
			case REDIS_REPLY_STRING:
				npn = tt(string(reply->str, reply->len));
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"GET %s:phase", to_string(f).c_str());
                if (reply == NULL) {
                    throw runtime_error("Error connecting to server");
                }
				phase = tt(string(reply->str, reply->len));
				freeReplyObject(reply);
				reply = (redisReply*)redisCommand(_rcontext, 
						"LRANGE %s:perm 0 -1", to_string(f).c_str());
                if (reply == NULL || reply->type != REDIS_REPLY_ARRAY) {
                    throw runtime_error("Error connecting to server");
                }
				for (auto i = 0u; i < reply->elements; i++) {
					auto preply = reply->element[i];
					auto permstring = string(preply->str, preply->len);
					perm.push_back(stoi(permstring));
				}
                if (perm.size() != tt_num_vars(f)) {
                    throw runtime_error("Error: perm.size() != tt_num_vars(f)");
                }
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				throw runtime_error(errorstring);
				break;
		}
		freeReplyObject(reply);
#endif
		return npn;
	}

}
