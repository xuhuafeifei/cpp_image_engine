add_rules("mode.debug", "mode.release")

add_requires("raylib")
add_requires("glm")

target("rmRenderer")
    set_kind("binary")
    add_files("src/*.cpp")
    add_packages("raylib")
    add_packages("glm")
--     add_packages("stb_image")
--     add_packages("thirdparty")
--
--     -- 添加 stb_image 头文件路径
--     add_includedirs("thirdparty/stb")

-- add_rules("mode.debug", "mode.release")
--
-- target("rmRenderer")
--     set_kind("binary")
--     add_files("src/*.cpp")
--     add_packages("raylib", {localdir = "E:\\c++_code\\rmRenderer\\lib\\raylib-5.5"})
--     add_packages("glm", {localdir = "E:\\c++_code\\rmRenderer\\lib\\glm"})
--
--     -- 指定 raylib 和 glm 的包含路径
--     add_includedirs("E:\\c++_code\\rmRenderer\\lib\\raylib-5.5\\src")
--     add_includedirs("E:\\c++_code\\rmRenderer\\lib\\")
--
--    -- 指定链接库
--    add_links("raylib")
--    add_links("glm")


-- xmake f -p windows -a x86_64 -m debug --toolchain=mingw

--
-- If you want to known more usage about xmake, please see https://xmake.io
--
-- ## FAQ
--
-- You can enter the project directory firstly before building project.
--
--   $ cd projectdir
--
-- 1. How to build project?
--
--   $ xmake
--
-- 2. How to configure project?
--
--   $ xmake f -p [macosx|linux|iphoneos ..] -a [x86_64|i386|arm64 ..] -m [debug|release]
--
-- 3. Where is the build output directory?
--
--   The default output directory is `./build` and you can configure the output directory.
--
--   $ xmake f -o outputdir
--   $ xmake
--
-- 4. How to run and debug target after building project?
--
--   $ xmake run [targetname]
--   $ xmake run -d [targetname]
--
-- 5. How to install target to the system directory or other output directory?
--
--   $ xmake install
--   $ xmake install -o installdir
--
-- 6. Add some frequently-used compilation flags in xmake.lua
--
-- @code
--    -- add debug and release modes
--    add_rules("mode.debug", "mode.release")
--
--    -- add macro definition
--    add_defines("NDEBUG", "_GNU_SOURCE=1")
--
--    -- set warning all as error
--    set_warnings("all", "error")
--
--    -- set language: c99, c++11
--    set_languages("c99", "c++11")
--
--    -- set optimization: none, faster, fastest, smallest
--    set_optimize("fastest")
--
--    -- add include search directories
--    add_includedirs("/usr/include", "/usr/local/include")
--
--    -- add link libraries and search directories
--    add_links("tbox")
--    add_linkdirs("/usr/local/lib", "/usr/lib")
--
--    -- add system link libraries
--    add_syslinks("z", "pthread")
--
--    -- add compilation and link flags
--    add_cxflags("-stdnolib", "-fno-strict-aliasing")
--    add_ldflags("-L/usr/local/lib", "-lpthread", {force = true})
--
-- @endcode
--
