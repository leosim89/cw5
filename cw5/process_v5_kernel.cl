uint vmin3(uint a, uint b, uint c)
{
	return min(a, min(b, c));
}

uint vmin4(uint a, uint b, uint c, uint d)
{
	return min(min(a, d), min(b, c));
}

uint vmin5(uint a, uint b, uint c, uint d, uint e)
{
	return min(e, min(min(a, d), min(b, c)));
}

void erode_pixel(
	unsigned l,
	unsigned x,
	unsigned levels,
	unsigned h,
	__global const uint *input,
	__global const uint *right,
	__global const uint *left,
	__global uint *output) {

	int index = l*(h + 1) + x;

	if (input[(l + 1) * h + l]){

		if (!left[(l + 1) * h + l]) {
			if (x == 0) 
				output[index + (h+1)] = vmin3(input[index], input[index + 1], right[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin3(input[index], input[index - 1], right[index]);
				output[index + (h + 2)] = 1;
			}
			else 
				output[index + (h + 1)] = vmin4(input[index], input[index + 1], input[index - 1], right[index]);
		}

		else if (!right[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmin3(input[index], input[index + 1], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin3(input[index], input[index - 1], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmin4(input[index], input[index + 1], input[index - 1], left[index]);
		}

		else {
			if (x == 0)
				output[index + (h + 1)] = vmin4(input[index], input[index + 1], right[index], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmin4(input[index], input[index - 1], right[index], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmin5(input[index], input[index + 1], input[index - 1], right[index], left[index]);
		}

	} 
	else {
		output[index + (h + 2)] = 0;
	}


}

uint vmax3(uint a, uint b, uint c)
{
	return max(a, max(b, c));
}

uint vmax4(uint a, uint b, uint c, uint d)
{
	return max(max(a, d), max(b, c));
}

uint vmax5(uint a, uint b, uint c, uint d, uint e)
{
	return max(e, max(max(a, d), max(b, c)));
}

void dilate_pixel(
	unsigned l,
	unsigned x,
	unsigned levels,
	unsigned h,
	__global const uint *input,
	__global const uint *right,
	__global const uint *left,
	__global uint *output) {

	int index = l*(h + 1) + x;

	if (input[(l + 1) * h + l]){

		if (!left[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmax3(input[index], input[index + 1], right[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax3(input[index], input[index - 1], right[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax4(input[index], input[index + 1], input[index - 1], right[index]);
		}

		else if (!right[(l + 1) * h + l]) {
			if (x == 0)
				output[index + (h + 1)] = vmax3(input[index], input[index + 1], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax3(input[index], input[index - 1], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax4(input[index], input[index + 1], input[index - 1], left[index]);
		}

		else {
			if (x == 0)
				output[index + (h + 1)] = vmax4(input[index], input[index + 1], right[index], left[index]);
			else if (x == h - 1){
				output[index + (h + 1)] = vmax4(input[index], input[index - 1], right[index], left[index]);
				output[index + (h + 2)] = 1;
			}
			else
				output[index + (h + 1)] = vmax5(input[index], input[index + 1], input[index - 1], right[index], left[index]);
		}

	}
	else {
		output[index + (h + 2)] = 0;
	}

}

__kernel void kernel_lx (
	int levels,
	unsigned h,
	__global const uint *input,
	__global const uint *right,
	__global const uint *left,
	__global uint *output)
{
	uint l=get_global_id(0);
	uint x=get_global_id(1);
	
	unsigned absl = abs(levels);
	if (levels<0) {
		erode_pixel(l, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
		dilate_pixel(l + absl, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
	} else {
		dilate_pixel(l, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
		erode_pixel(l + absl, x, absl, h, &input[0], &right[0], &left[0], &output[0]);
	}
}