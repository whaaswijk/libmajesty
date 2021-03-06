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

	function_store::function_store(const std::string& server_url, const unsigned port) {
		auto path = boost::filesystem::unique_path();
		_logfile = path.native();
#ifndef _WIN32
		_rcontext = redisConnect(server_url.c_str(), port);
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
		return min_size_xmg(f, 0).get();
	}


	optional<unsigned> function_store::get_last_size(const cirkit::tt& f) {
		optional<unsigned> last_size;
#ifndef _WIN32
		redisReply* reply = (redisReply*)redisCommand(_rcontext, "GET %s:last_size", to_string(f).c_str());
		if (reply == NULL) {
			throw runtime_error("Error connecting to server");
		}
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				// Haven't timed out on this function
				break;
			case REDIS_REPLY_STRING:
				last_size = unsigned(stoi(string(reply->str, reply->len)));
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
		}
#endif
		return last_size;
	}

	void function_store::set_heuristic_size_xmg(const cirkit::tt& f, const string& expr, unsigned timeout) {
#ifndef _WIN32
		redisReply* reply = (redisReply*)redisCommand(_rcontext, "SET %s:heuristic %s",
				to_string(f).c_str(), expr.c_str());
		if (reply == NULL) {
			throw runtime_error("Error connecting to server");
		}
		freeReplyObject(reply);
		reply = (redisReply*)redisCommand(_rcontext, "SET %s:heuristic_timeout %u",
				to_string(f).c_str(), timeout);
		if (reply == NULL) {
			throw runtime_error("Error connecting to server");
		}
		freeReplyObject(reply);
#endif
	}
			
	optional<string> function_store::heuristic_size_xmg(const cirkit::tt& f, unsigned timeout) {
#ifndef _WIN32
		// Make sure that we've timed out for this function before. If not, no heurstic version exists.
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "GET %s:heuristic_timeout", to_string(f).c_str());
		if (reply->type == REDIS_REPLY_STRING) {
			// We only attempt to synthesize again if our timeout is bigger this time
			auto oldtimeout = unsigned(stoi(string(reply->str, reply->len)));
			freeReplyObject(reply);
			if (oldtimeout < timeout) {
				// If our timeout is now bigger we want to see if we can get a
				// better heurstic version with this new timeout
				return boost::none;
			}
			// Return the current heuristic
			reply = (redisReply*) redisCommand(_rcontext, "GET %s:heuristic", to_string(f).c_str());
			if (reply->type == REDIS_REPLY_STRING) {
				auto expr = string(reply->str, reply->len);
				freeReplyObject(reply);
				return expr;
			} else if (reply->type == REDIS_REPLY_NIL) {
				// No heuristic was computed yet
				freeReplyObject(reply);
				return boost::none;
			} else {
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
			}
		} else if (reply->type == REDIS_REPLY_NIL) {
			freeReplyObject(reply);
			return boost::none;
		} else {
			auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
			auto errorstring = errorformat.str();
			freeReplyObject(reply);
			throw runtime_error(errorstring);
		}
#else
		return "placeholder";
#endif
	}

	optional<string> function_store::min_size_xmg(const cirkit::tt& f, unsigned timeout) {
		// First see if the optimum xmg has already been computed
		optional<string> expr;
#ifndef _WIN32
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "GET %s:expr", to_string(f).c_str());
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				// Maybe the expression doesn't exist because we've timed out trying to compute it
				reply = (redisReply*) redisCommand(_rcontext, "GET %s:timeout", to_string(f).c_str());
				if (reply->type == REDIS_REPLY_STRING) {
					// We only attempt to synthesize again if our timeout is bigger this time
					auto oldtimeout = unsigned(stoi(string(reply->str, reply->len)));
					freeReplyObject(reply);
					if (oldtimeout >= timeout) {
						cout << "Skipping because of previous timeout" << endl;
						break;
					}
				} else {
					assert(reply->type == REDIS_REPLY_NIL);
					freeReplyObject(reply);
				}
				expr = exact_xmg_expression(f, timeout);
				if (expr) {
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:expr %s",
							to_string(f).c_str(), expr.get().c_str());
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
				} else { // Ensure we don't try to synthesize it again with the same timeout
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:timeout %u",
						to_string(f).c_str(), timeout);
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
					auto olast_size = last_size_from_file("cirkit.log");
					assert(olast_size);
					auto last_size = olast_size.get();
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:last_size %u",
						to_string(f).c_str(), last_size);
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
				}
				break;
			case REDIS_REPLY_STRING:
				expr = string(reply->str, reply->len);
				freeReplyObject(reply);
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
				break;
		}
#else
		expr = "<ab0>";
#endif
		return expr;
	}
	
	tt function_store::npn_canon(const tt& f, tt& phase, vector<unsigned>& perm) {
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
	
	/*
	boost::optional<std::string> function_store::get_size_optimum_ntk_ns(const tt& f, synth_spec* spec, unsigned conflict_limit) {
		// First see if the optimum xmg has already been computed
		optional<string> expr;
#ifndef _WIN32
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "GET %s:expr", to_string(f).c_str());
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				freeReplyObject(reply);
				// Maybe the expression doesn't exist because we've reached the
				// conflict limit trying to compute it 
				reply = (redisReply*) redisCommand(_rcontext, "GET %s:limit", to_string(f).c_str());
				if (reply->type == REDIS_REPLY_STRING) {
					// We only attempt to synthesize again if our conflict limit is bigger this time
					auto oldlimit = unsigned(stoi(string(reply->str, reply->len)));
					freeReplyObject(reply);
					if (oldlimit >= conflict_limit) {
						cout << "Skipping because of previous limit breach" << endl;
						break;
					}
				} else {
					assert(reply->type == REDIS_REPLY_NIL);
					freeReplyObject(reply);
				}
				//auto opt_size_ntk = 
				expr = logic_ntk_to_string(size_optimum_ntk_ns<CMSat::SATSolver>(f, spec));
				if (expr) {
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:expr %s",
							to_string(f).c_str(), expr.get().c_str());
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
				} else { // Ensure we don't try to synthesize it again with the same conflict limit
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:limit %u",
						to_string(f).c_str(), conflict_limit);
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
					auto olast_size = last_size_from_file("cirkit.log");
					assert(olast_size);
					auto last_size = olast_size.get();
					reply = (redisReply*)redisCommand(_rcontext, "SET %s:last_size %u",
						to_string(f).c_str(), last_size);
					if (reply == NULL) {
						throw runtime_error("Error connecting to server");
					}
					freeReplyObject(reply);
				}
				break;
			case REDIS_REPLY_STRING:
				expr = string(reply->str, reply->len);
				freeReplyObject(reply);
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
				break;
		}
#else
		expr = "<ab0>";
#endif
		return expr;

	}
	*/

	store_stats function_store::get_store_stats() {
		store_stats stats;
		
#ifndef _WIN32
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "KEYS *");
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL:
				cout << "No keys found" << endl;
				break;
			case REDIS_REPLY_ARRAY:
				cout << "Nr keys: " << reply->elements << endl;
				for (auto i = 0u; i < reply->elements; i++) {
					auto element = reply->element[i];
					//auto key = string(reply->str, reply->len);
					cout << "reply: " << element->type << endl;
				}
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
				break;
		}
#endif
		return stats;
	}
			
	void function_store::set_entry(const cirkit::tt& function, const std::string& expr, const bool is_exact, const unsigned conflict_limit) {
		function_store_entry entry(expr, is_exact, conflict_limit);
		set_entry(function, entry);
	}

	void function_store::set_entry(const cirkit::tt& function, const function_store_entry& entry) {
#ifndef _WIN32
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "SET %s %s", to_string(function).c_str(), entry.to_string().c_str());
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		freeReplyObject(reply);
#else
#endif
	}

	boost::optional<function_store_entry> function_store::get_entry(const cirkit::tt& function) {
		function_store_entry entry;
		string entry_string;
#ifndef _WIN32
		redisReply* reply = (redisReply*) redisCommand(_rcontext, "GET %s", to_string(function).c_str());
        if (reply == NULL) {
            throw runtime_error("Error connecting to server");
        }
		switch (reply->type) {
			case REDIS_REPLY_NIL: // This entry doesn't exist yet
				freeReplyObject(reply);
				return boost::none;
				break;
			case REDIS_REPLY_STRING:
				entry_string = string(reply->str, reply->len);
				entry = to_function_store_entry(entry_string);
				freeReplyObject(reply);
				return std::move(entry);
				break;
			default:
				auto errorformat = boost::format("Unable to handle reply type: %s") % reply->type;
				auto errorstring = errorformat.str();
				freeReplyObject(reply);
				throw runtime_error(errorstring);
				break;
		}
#endif
		return std::move(entry);
	}
	
	// Function store entry format: expr-is_exact-conflict_limit
	function_store_entry to_function_store_entry(const std::string& str) {
		function_store_entry entry;

		const auto tokens = split(str, '-');
		const auto ntokens = tokens.size();
		assert(ntokens == 3);

		entry.expr = tokens[0];
		entry.is_exact = stoi(tokens[1]);
		entry.conflict_limit = stoi(tokens[2]);

		return entry;
	}
}
