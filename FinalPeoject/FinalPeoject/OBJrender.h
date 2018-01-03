#pragma once

#ifndef __OBJrender_H__
#define __OBJrender_H__

#include <array>
#include <functional>
#include <memory>
#include <vector>
#include "VBOobj.h"

#define any_layer -1
#define world_layer 0
#define ui_layer 1
#define world_transparent_layer 2


//typedef std::array<float, 3> mycolor3;
typedef std::array<float, 4> mycolor4;

namespace OBJrender
{
	template <typename T>
	using is_vector = std::is_same<T, std::vector<typename T::value_type, typename T::allocator_type>>;

	template <typename T, typename = void>
	struct dodelete
	{
		void operator()(T& obj) const
		{
		}
	};
	template <typename T>
	struct dodelete<T, std::enable_if_t<std::is_pointer<T>::value>>
	{
		void operator()(T& obj) const
		{
			delete obj;
		}
	};
	template <typename T>
	struct dodelete<T, std::enable_if_t<is_vector<T>::value>>
	{
		void operator()(T& obj) const
		{
			dodelete<T::value_type> deleter;
			for (auto&& i : obj)
			{
				deleter(i);
			}
		}
	};

	template <typename T>
	class OBJManger
	{
	public:
		typedef T type;
		typedef T& rtype;
		typedef const T& crtype;
		std::shared_ptr<T> obj;

		OBJManger() = delete;

		virtual ~OBJManger()
		{
			if (obj.use_count() == 1)
			{
				dodelete<T>()(*obj);
			}			
		};

		OBJManger(const std::shared_ptr<T>& O) : obj(O){}
	};
//////////////////////////////////////////////////////////////////////////
	template <typename T>
	struct settingVBO
	{
		void operator()(T& obj)
		{
			setting(obj);
		}
		void setting(T& obj);
	};
	
	template <typename T>
	struct displayfunc
	{
		void operator()(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const
		{
			display(obj, colorf, il);
		}

		void display(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const;
	};

	/*template <typename T>
	struct uidisplayfunc
	{
		void operator()(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const
		{
			display(obj, colorf, il);
		}

		void display(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const;
	};*/
//////////////////////////////////////////////////////////////////////////
	/*class Render
	{
	public:
		Render(const mycolor4& colorf_ = {0, 0, 0, 1}) : colorf(colorf_){}
		virtual ~Render(){};

		virtual void Display() const = 0;

		void setColorf(float R_, float G_, float B_, float A_ = 1.0f)
		{
			colorf = {R_, G_, B_, A_};
		};

	protected:
		mycolor4 colorf;
	};

	template <typename T, typename dfn = displayfunc<T>>
	class Normalrender :public Render, public OBJManger<T>
	{
	public:
		template <typename T2>
		Normalrender(T2&& O, const mycolor4& colorf_ = {0, 0, 0, 1}) : Normalrender(std::make_shared<T>(std::forward<T2>(O)), colorf_)
		{
		}
		Normalrender(std::shared_ptr<T> O, const mycolor4& colorf_ = {0, 0, 0, 1}) : Render(colorf_),  OBJManger(O)
		{
		}

		virtual void Display() const
		{
			dis(*obj, colorf);
		}		
		template<typename Fn>
		void Display_fn(Fn fn) const
		{
			fn(*obj, colorf);
		}
	protected:
		dfn dis;
	};

	template <typename T, typename dfn = displayfunc<T>, typename sfn = settingVBO<T>>
	class VBOrender :public Normalrender<T, dfn>
	{
	public:
		template<typename T2>
		VBOrender(T2&& O, const mycolor4& colorf_ = {0, 0, 0, 1}) : VBOrender(std::make_shared<T>(std::forward<T2>(O)), colorf_)
		{
		}
		VBOrender(std::shared_ptr<T> O, const mycolor4& colorf_ = {0, 0, 0, 1}) : Normalrender(O, colorf_)
		{
			genVBO(obj.get());
			sett(*obj);
		}

		virtual void Display() const
		{
			dis(*obj, colorf);
		}
		template<typename Fn>
		void Display_fn(Fn fn) const
		{
			fn(*obj, colorf);
		}
	protected:
		sfn sett;
	};

	template<typename T, typename dfn = displayfunc<T>, typename sfn = settingVBO<T>>
	using template_Render = std::conditional_t<std::is_base_of<VBOobj, T>::value, VBOrender<T, dfn, sfn>, Normalrender<T, dfn>>;

	class uiRender_base
	{
	public:
		uiRender_base(const mycolor4& uicolorf_ = {0, 0, 0, 1}) : uicolorf(uicolorf_){}
		virtual void uiDisplay() const = 0;

		void setuiColorf(float R_, float G_, float B_, float A_ = 1.0f)
		{
			uicolorf = {R_, G_, B_, A_};
		};
	protected:
		mycolor4 uicolorf;
	};

	template<typename T, typename udfn = uidisplayfunc<T>>
	class uiRender :public template_Render<T>, public uiRender_base
	{
	public:
		uiRender(T*&& O, const mycolor4& colorf_ = {0, 0, 0, 1}, const mycolor4& uicolorf_ = {0, 0, 0, 1}) : template_Render<T>(std::forward<T*>(O), colorf_), uiRender_base(uicolorf_)
		{
		}
		uiRender(const std::shared_ptr<T>& O, const mycolor4& colorf_ = {0, 0, 0, 1}, const mycolor4& uicolorf_ = {0, 0, 0, 1}) : template_Render<T>(O, colorf_), uiRender_base(uicolorf_)
		{
		}
		virtual void Display() const
		{
			dis(*obj, colorf);
		}
		virtual void uiDisplay() const
		{
			uidis(*obj, uicolorf);
		}
	protected:
		udfn uidis;
	};*/
	//////////////////////////////////////////////////////////////////////////

	template<typename T, typename = void>
	struct setter
	{
		setter(std::shared_ptr<T> obj)
		{
		}
		bool _releaseVBO(std::shared_ptr<T> obj)
		{
			return false;
		}
	};

	template<typename T>
	struct setter<T, std::enable_if_t<std::is_base_of<VBOobj, T>::value>>
	{
		typedef settingVBO<T> sfn;
		setter(std::shared_ptr<T> obj)
		{
			//std::cout << "set " << __FUNCDNAME__ << std::endl;
			genVBO(obj.get());
			sett(*obj);
		}
		sfn sett;

		bool _releaseVBO(std::shared_ptr<T> obj)
		{
			return obj->releaseVBO();
		}
	};

	template<typename T, int layer>
	struct displayfunc2
	{
		void operator()(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const
		{
			display(obj, colorf, il);
		}

		void display(const T& obj, const mycolor4& colorf = {0, 0, 0, 1}, std::initializer_list<float> il = {}) const;
	};

	template<typename T, int layer>
	using displayfunc3 = std::conditional_t<layer == any_layer, displayfunc<T>, displayfunc2<T, layer>>;

	/*namespace test
	{
		template<int l>
		class layers;

		template<>
		class layers<wolrd_layer> :public Render
		{
		public:
			using Render::colorf;
			layers<wolrd_layer>(const mycolor4& c = {0,0,0,1}) : Render(c)
			{
			}
		};

		template<>
		class layers<ui_layer> :public uiRender_base
		{
		public:
			layers<ui_layer>(const mycolor4& c = {0,0,0,1}) : uiRender_base(c), colorf(c)
			{
			}
			virtual void uiDisplay() const
			{}
			mycolor4 colorf;
		};

		template<typename T, int... layer>
		class test_render :public layers<layer>..., public OBJManger<T>, public setter<T>
		{
		public:
			template<typename T2, typename... T3>
			test_render(T2&& O, const T3&... colorf_) : test_render(std::make_shared<T>(std::forward<T2>(O)), colorf_...)
			{
			}
			template<typename... T3>
			test_render(std::shared_ptr<T> O, const T3&... colorf_) : layers<layer>(colorf_)..., OBJManger(O), setter(O), colorfarr({colorf_...})
			{
			}
			template<typename T2>
			test_render(T2&& O) : test_render(std::make_shared<T>(std::forward<T2>(O)))
			{
			}
			test_render(std::shared_ptr<T> O) : layers<layer>()..., OBJManger(O), setter(O), colorfarr({{0,0,0,1}})
			{
			}

			virtual void Display() const
			{
				std::get<0>(dis)(*obj, colorfarr[0]);
			}
			template<int l>
			void Display_l() const
			{
				std::get<displayfunc2<T, l>>(dis)(*obj, layers<l>::colorf);
			}
			template<typename Fn>
			void Display_fn(Fn fn) const
			{
				fn(*obj, colorf);
			}
		protected:
			std::tuple<displayfunc2<T, layer>...> dis;
			std::array<mycolor4, sizeof...(layer)> colorfarr;
		};
	}*/

	namespace test2
	{
		typedef std::function<bool()> bfun;

		static bfun defaultbun = []()
		{
			return true;
		};

		class clickableui
		{
		public:
			bfun checkfun = defaultbun;
			bool operator()() const
			{
				return checkfun();
			}
		};

		class moveableui
		{
		public:
			double x;
			double y;
		};

		template<int l>
		class layerbase
		{
		public:
			virtual ~layerbase() {};
			virtual void Display() const = 0;
			virtual bool checkclick() const = 0;
			virtual bool releaseVBO() = 0;
		};

		template<typename T, int l = any_layer>
		class layerrender :virtual public OBJManger<T>, virtual public setter<T>, public layerbase<l>
		{
		public:
			template<typename T2>
			layerrender(T2&& O, const mycolor4& colorf_ = {0,0,0,1}) : layerrender(std::make_shared<T>(std::forward<T2>(O)), colorf_)
			{
			}
			layerrender(std::shared_ptr<T> O, const mycolor4& colorf_ = {0,0,0,1}) : OBJManger(O), setter(O), colorf(colorf_), checkfun(defaultbun)
			{
			}
			
			virtual void Display() const final
			{
				//std::cout << __FUNCDNAME__ << std::endl;
				//std::cout << "layer " << l << " obj.count: " << obj.use_count() <<  std::endl;
				if (checkfun())
				{
					dis(*obj, colorf);
				}
			}
			
			void setbfun(bfun b)
			{
				checkfun = b;
			}

			virtual bool releaseVBO()
			{
				return _releaseVBO(obj);
			}

			template<typename tt, typename enable = void>
			struct ctttt
			{
				bool operator()(T* obj)
				{
					return false;
				}
			};

			template<typename tt>
			struct ctttt<tt, std::enable_if_t<std::is_base_of<clickableui, tt>::value>>
			{
				bool operator()(T* obj)
				{
					return obj->checkfun();
				}
			};

			virtual bool checkclick() const final
			{
// 				if constexpr(std::is_base_of<clickable, T>::value)
// 				{
// 					return obj->checkfun();
// 				}
// 				return false;

				return ctttt<T>()(obj.get());
			}

			template<typename Fn>
			void Display_fn(Fn fn, const mycolor4& colorf) const
			{
				fn(*obj, colorf);
			}
		protected:
			displayfunc3<T, l> dis;
			mycolor4 colorf;
			bfun checkfun;
		};

		template<typename T, int... layer>
		class multilayerrender :public layerrender<T, layer>...
		{
		public:
			template<typename T2, typename... T3>
			multilayerrender(T2&& O, const T3&... colorf_) : multilayerrender(std::make_shared<T>(std::forward<T2>(O)), colorf_...)
			{
			}
			template<typename... T3>
			multilayerrender(std::shared_ptr<T> O, const T3&... colorf_) :OBJManger(O), setter(O), layerrender<T, layer>(O, colorf_)...
			{
			}
			template<typename T2>
			multilayerrender(T2&& O) : multilayerrender(std::make_shared<T>(std::forward<T2>(O)))
			{
			}
			multilayerrender(std::shared_ptr<T> O) :OBJManger(O), setter(O), layerrender<T, layer>(O)...
			{
			}

			template<int l>
			void Display_l() const
			{
				layerrender<T, l>::Display();
			}
			template<int l>
			void setbfun_l(bfun b)
			{
				layerrender<T, l>::checkfun = b;
			}
			template<typename Fn>
			void Display_fn(Fn fn, const mycolor4& colorf) const
			{
				fn(*obj, colorf);
			}

			virtual bool releaseVBO()
			{
				return _releaseVBO(obj);
			}
		};
	}
}

#define _mydisplay(A) OBJrender::displayfunc<OBJrender::A::type>::display(A::crtype obj, const mycolor4& colorf, std::initializer_list<float> il) const
#define mysetter(A) OBJrender::settingVBO<OBJrender::A::type>::setting(A::rtype obj)

#define mydisplay2(A, layer) OBJrender::displayfunc2<OBJrender::A::type, layer>::display(A::crtype obj, const mycolor4& colorf, std::initializer_list<float> il) const
#define mydisplay(A) mydisplay2(A, world_layer)
#define myuidisplay(A) mydisplay2(A, ui_layer)
#define mydisplay_t(A) mydisplay2(A, world_transparent_layer)

#endif