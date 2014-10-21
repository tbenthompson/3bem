#ifndef __FMM_OPENCL_H
#define __FMM_OPENCL_H
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
 
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

std::string get_platform_name(cl_platform_id id) {
    size_t size = 0;
    clGetPlatformInfo (id, CL_PLATFORM_NAME, 0, nullptr, &size);
    std::string result;
    result.resize(size);
    clGetPlatformInfo(id, CL_PLATFORM_NAME, size,
                      const_cast<char*> (result.data ()),
                      nullptr);
    return result;
}

std::string get_device_name (cl_device_id id) {
    size_t size = 0;
    clGetDeviceInfo(id, CL_DEVICE_NAME, 0, nullptr, &size);
    std::string result;
    result.resize(size);
    clGetDeviceInfo(id, CL_DEVICE_NAME, size,
                    const_cast<char*> (result.data ()),
                    nullptr);
    return result;
}

void check_error(cl_int error)
{
    if (error != CL_SUCCESS) {
        std::cerr << "OpenCL call failed with error " << error << std::endl;
        std::exit(1);
    }
}

std::string load_kernel(const char* name) {
    std::ifstream file(name);
    std::string result((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    return result;
}

cl_program create_program(const std::string& source, cl_context context)
{
    size_t lengths[1] = {source.size()};
    const char* sources[1] = {source.data()};
    cl_int error = 0;
    cl_program program = clCreateProgramWithSource(context, 1, 
                                    sources, lengths, &error);
    check_error(error);
    return program;
}

void setup_ocl() {
    cl_uint platform_count = 0;
    clGetPlatformIDs(0, nullptr, &platform_count);

    if (platform_count == 0) {
        std::cerr << "No OpenCL platform found" << std::endl;
        std::exit(1);
    } else {
        std::cout << "Found " << platform_count << " platform(s)" << std::endl;
    }

    std::vector<cl_platform_id> platform_ids(platform_count);
    clGetPlatformIDs(platform_count, platform_ids.data(), nullptr);

    for (cl_uint i = 0; i < platform_count; ++i) {
        std::cout << "\t (" << (i+1) << ") : " << 
            get_platform_name(platform_ids[i]) << std::endl;
    }

    int platform_choice = 0;

    // http://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/clGetDeviceIDs.html
    cl_uint device_count = 0;
    clGetDeviceIDs(platform_ids[platform_choice], CL_DEVICE_TYPE_ALL, 
                   0, nullptr, &device_count);

    if (device_count == 0) {
        std::cerr << "No OpenCL devices found" << std::endl;
        std::exit(1);
    } else {
        std::cout << "Found " << device_count << " device(s)" << std::endl;
    }

    std::vector<cl_device_id> device_ids(device_count);
    clGetDeviceIDs(platform_ids[platform_choice],
                   CL_DEVICE_TYPE_ALL, device_count,
                   device_ids.data(), nullptr);

    for (cl_uint i = 0; i < device_count; ++i) {
        std::cout << "\t (" << (i+1) << ") : " << 
                 get_device_name(device_ids[i]) << std::endl;
    }

    const cl_context_properties contextProperties [] =
    {
            CL_CONTEXT_PLATFORM,
            reinterpret_cast<cl_context_properties>(platform_ids[platform_choice]),
            0, 0
    };

    /* Create OpenCL context */
    cl_int error = CL_SUCCESS; 
    cl_context context = clCreateContext(contextProperties, device_count,
		                          device_ids.data(), nullptr, nullptr, &error);
    check_error(error);
     
    /* Create Command Queue */
    cl_command_queue command_queue = clCreateCommandQueue(context, device_ids[0], 0, &error);
    check_error(error);
     
    std::string source = load_kernel("src/kernels/abc.cl");
    std::cout << source << std::endl;
    cl_program program = create_program(source, context);

    unsigned int n = 1000;
    cl_mem memobj = clCreateBuffer(context, CL_MEM_READ_WRITE, n * sizeof(float), nullptr, &error);
    check_error(error);

    error = clBuildProgram(program, device_count, device_ids.data(), "", nullptr, nullptr);
    check_error(error);

    cl_kernel kernel = clCreateKernel(program, "vecAdd", &error);
    check_error(error);

    error = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&memobj);
    check_error(error);

    // error = clEnqueueTask(command_queue, kernel, 0, nullptr, nullptr);
    // check_error(error);
    size_t global_size[1] = {n};
    error = clEnqueueNDRangeKernel(command_queue, kernel, 1, nullptr,
                                   global_size, nullptr, 0, nullptr, nullptr);

    float* result = new float[n];
    error = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0, n * sizeof(float),
                                result, 0, nullptr, nullptr); 
    for (unsigned int i = 0; i < n; i++) {
        // std::cout << result[i] << std::endl;
    }

    error = clFlush(command_queue);
    check_error(error);
    error = clFinish(command_queue);
    check_error(error);
    error = clReleaseKernel(kernel);
    check_error(error);
    error = clReleaseProgram(program);
    check_error(error);
    error = clReleaseMemObject(memobj);
    check_error(error);
    error = clReleaseCommandQueue(command_queue);
    check_error(error);
    error = clReleaseContext(context);
    check_error(error);
// 
//     cl_mem memobj = NULL;
//     cl_program program = NULL;
//     cl_kernel kernel = NULL;
// 
//     cl_platform_id platform_id = NULL;;
//         /* Create Memory Buffer */
//     memobj = clCreateBuffer(context, CL_MEM_READ_WRITE,MEM_SIZE * sizeof(char), NULL, &ret);
//      
//     /* Create Kernel Program from the source */
//     std::string source = load_kernel("abc.cl");
//     program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
//     (const size_t *)&source_size, &ret);
//      
//     /* Build Kernel Program */
//     ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
//      
//     /* Create OpenCL Kernel */
//     kernel = clCreateKernel(program, "hello", &ret);
//      
//     /* Set OpenCL Kernel Parameters */
//     ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);
//      
//     /* Execute OpenCL Kernel */
//     ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);
//      
//     /* Copy results from the memory buffer */
//     ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
//     MEM_SIZE * sizeof(char),string, 0, NULL, NULL);
//      
//     /* Display Result */
//     puts(string);
//      
//     /* Finalization */
//     ret = clFlush(command_queue);
//     ret = clFinish(command_queue);
//     ret = clReleaseKernel(kernel);
//     ret = clReleaseProgram(program);
//     ret = clReleaseMemObject(memobj);
//     ret = clReleaseCommandQueue(command_queue);
//     ret = clReleaseContext(context);
}
#endif
