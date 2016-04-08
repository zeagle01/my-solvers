template<class T>
class ExpansionFunction{
public:
	virtual T operator()(T x, T y) = 0;
};


template<class T> class UniformExpansion :public ExpansionFunction<T>{
public:
	virtual T operator()(T x, T y){
		return x + y;
	}
};

template<class T> class print{

public:
	void operator()(T e){
		cout << e << endl;
	}
};


