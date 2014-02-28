all:
	g++ tools/MergeFusion/main.cpp tools/MergeFusion/Merger.cpp tools/MergeFusion/StandardFusion.cpp tools/MergeFusion/common.cpp -o MergeFusion

clean:
	rm MergeFusion
