#include <algorithm>

class ArgumentParser{
    private:
        std::vector <std::string> tokens;
    public:
        ArgumentParser(int argc, char* argv[]){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        bool findOption(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option) != this->tokens.end();
        }

        void getArgument(const std::string &option, std::string& arg) const{
            static const std::string empty_string("");
            arg = empty_string;
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                arg = *itr;
            }
        }

        void getArgument(const std::string &option, int& arg) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                arg = stoi(*itr);
            }
        }
        
        void getArgument(const std::string &option, std::vector<std::string> &args) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            for (itr != this->tokens.end(); itr != this->tokens.end(); ++itr){
                string s = *itr;
                if(s==option) continue;
                if(s.find("-")==0) break;
                args.push_back(s);
            }
        }
};