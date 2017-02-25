#include <string>
#include <boost/optional.hpp>

namespace majesty
{
	boost::optional<std::string>
		hex_to_bool_string(const std::string& hex_str)
	{
		std::string res;

		for (auto i = 2u; i < hex_str.size(); i++)
		{
			const auto& hex_char = hex_str[i];
			char intval;
			if (hex_char >= '0' && hex_char <= '9')
			{
				intval = hex_char - '0';
			}
			else if (hex_char >= 'a' && hex_char <= 'f')
			{
				intval = (hex_char - 'a') + 10;
			}
			else if (hex_char >= 'A' && hex_char <= 'F')
			{
				intval = (hex_char - 'A') + 10;
			}
			else
			{
				return boost::none;
			}
			for (auto i = 3; i >= 0; i--)
			{
				res.push_back( ((intval >> i) & 1) == 0 ? '0' : '1' );
			}
		}
		
		return res;
	}
	
}

