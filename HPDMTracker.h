/*
 *  HPDMTracker.h
 *  CurveMatching
 *
 *  Created by roy_shilkrot on 2/15/13.
 *  Copyright (c) 2013 MIT
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * from "Wormholes in Shape Space: Tracking through Discontinuous Changes in Shape", Heap & Hogg, 1998
 */

#include "AbstractAlgorithm.h"

#include "HPDM.h"
#include "MathGLTools.h"

class HPDMTracker : public AbstractAlgorithm {
public:
    class Particle {
	public:
		//according to section 3.1
		Point2f centroid;
		float scale;
		float orientation;
		Mat_<float> deformation;
		int patch_num;
		vector<pair<int,Point2f> > dx_dy;
        
        Particle() {}
        
        //cpy c'tor
        Particle(const Particle& other):
        centroid(other.centroid),
        scale(other.scale),
        orientation(other.orientation),
        patch_num(other.patch_num)
        {
            other.deformation.copyTo(deformation);
            dx_dy.insert(dx_dy.begin(), other.dx_dy.begin(), other.dx_dy.end());
        }
	};
private:
	
	HPDM<float> hpdm;
	vector<vector<Point2f> > m_curves;
	Mat_<double> TransMat;
	Mat_<double> Cumulative;
	vector<Particle> particles;
	int num_particles;
	vector<float> unnormalized_weights;
	vector<float> normalized_weights;
	vector<float> cumulative_sum;
	vector<Point2f> mean_shape;
    Particle mean_particle;
    
	int curve_length;
    cv::Size bounds;
    int m_iter;
    float initscale;
    float N_thr;
    float N_eff;
    
	RNG rng_;

	void buildTransitionMatrix() {
		Mat show(480,640,CV_8UC3);
		//--create transition matrix from sequential samples
		TransMat.create(hpdm.getNumPatches(),hpdm.getNumPatches());
		TransMat.setTo(0);
		for (int i=0; i<m_curves.size()-1; i++) {
			
			//get nearest patch for i and i-1
			Mat_<double> dists_i_1, dists_i; 
			int pi_1 = hpdm.nearestPatch(m_curves[i], dists_i_1);
			int pi = hpdm.nearestPatch(m_curves[i+1], dists_i);
			Mat_<double> Ti_i_1 = dists_i_1.t() * dists_i;
//			TransMat += Ti_i_1; //eqn. (1)
			
			//Use only nearest patch relationship
            if(pi_1 >= 0 && pi >= 0) {
                TransMat(pi_1,pi)++;
//			TransMat(pi,pi_1)++;
            }
			
			
//			 show.setTo(0);
//            vector<Point2f> curve1(m_curves[i].begin(),m_curves[i].end());
//            Mat(curve1) += Scalar(150,150);
//			 drawOpenCurve(show, curve1, Scalar(255), 1);
//            curve1.clear(); curve1.insert(curve1.begin(),m_curves[i+1].begin(),m_curves[i+1].end());
//            Mat(curve1) += Scalar(150,150);
//			 drawOpenCurve(show, curve1, Scalar(255), 1);
//
//            stringstream ss;
//			 ss << "patch " << pi_1 << " -> patch " << pi;
//            putText(show, ss.str(), cv::Point(10,10), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(0,0,255), 1);
//			 imshow("curve",show);
//			 
//			 double maxVal,minVal;
//			 minMaxLoc(Ti_i_1, &minVal, &maxVal);
//			 cout << "maxval " << maxVal << "\n";
//			 imshow_<276,276>("t mat", (Ti_i_1 - minVal) / (maxVal-minVal));
//			 
//			 waitKey(0);
            
			 
		}
		destroyAllWindows();
		
		//Normalize Transition Matrix to create a PDF and build Cumulative Distribution Function matrix
		Cumulative.create(hpdm.getNumPatches(),hpdm.getNumPatches());
		for (int i=0; i<TransMat.rows; i++) {
			Scalar sum_row = sum(TransMat.row(i));
			
			//eqn. (2)
			TransMat.row(i) /= sum_row[0];
			
			//eqn. (3)
			Cumulative(i,0) = TransMat(i,0);
			for (int j=1; j<TransMat.cols; j++) {
				Cumulative(i,j) = Cumulative(i,j-1) + TransMat(i,j);
			}
		}
		
        
        /*
		double maxVal,minVal;
		minMaxLoc(TransMat, &minVal, &maxVal);
		cout << "maxval " << maxVal << "\n";
		imshow_<276,276>("t mat", (TransMat - minVal) / (maxVal-minVal));
		minMaxLoc(Cumulative, &minVal, &maxVal);
		cout << "maxval " << maxVal << "\n";
		imshow_<276,276>("Cumulative", (Cumulative - minVal) / (maxVal-minVal));
		waitKey(0);
        */
        
	}
	
	void initializeParticles() {
		//randomize N particles
		particles.resize(num_particles);
		unnormalized_weights.resize(num_particles,1.0f/(float)num_particles);
		normalized_weights.resize(num_particles,1.0f/(float)num_particles);
		cumulative_sum.resize(num_particles,0);
        normalizeWeights();
        
		for (int i=0; i<num_particles; i++) {
			//randomize
			particles[i].orientation = rng_.uniform(80.0,100.0);
			particles[i].centroid.x = rng_.uniform(100,540);
			particles[i].centroid.y = rng_.uniform(100,380);
			particles[i].scale = rng_.uniform(1.0,1.5);
			particles[i].patch_num = rng_.uniform(0, hpdm.getNumPatches()-1);
			hpdm.getPatchMean(particles[i].patch_num).copyTo(particles[i].deformation);
//			particles[i].deformation.create(1,hpdm.getShapeSpaceDim());
//			rng_.fill(particles[i].deformation, RNG::NORMAL, 0, 0.33);
		}
	}
	
public:
	float calculateFitness(Particle& p, const Mat& dx, const Mat& skin_mask, vector<Point2f>& curve, int curve_step = 2, float step_scale = 1.0) {
//		//get rebuilt curve
//		vector<Point2f> curve = hpdm.rebuildCurve(p.deformation);
//		//transform according to particle
//		Mat_<float> Tpca = Mat_<float>::eye(3,3);
//		getRotationMatrix2D(Point2f(0,0), p.orientation, p.scale).copyTo(Tpca.rowRange(0,2));
//		Tpca = (Mat_<float>(3,3) << 1,0,p.centroid.x/p.scale, 0,1,p.centroid.y/p.scale, 0,0,1) * Tpca;
//		//	cout << "rot mat from angle\n"<<getRotationMatrix2D(cp, angle * 180.0 / CV_PI, 1.0) << "\n";
////		cout << "trans\n"<<Tpca<<"\n";
//		transform(curve,curve,Tpca.rowRange(0,2));
		int steps = max(10,cvRound(40 * step_scale));
		int low = -floor(steps/2),high = ceil(steps/2);
		steps = high-low;

		curve = getShapeForParticle(p);
		p.dx_dy.clear();
		vector<float> edge_history;
#undef _DEBUG_
#ifdef _DEBUG_
		vector<Mat> splt; split(dx, splt);
//		Mat gray8; splt[0].convertTo(gray8, CV_8UC1, 255);
//		Mat new_gray(gray8.size(),CV_8UC1,Scalar::all(0)); gray8.copyTo(new_gray);
		Mat show; //cvtColor(new_gray, show, CV_GRAY2BGR);
        splt[0] = Mat::zeros(splt[0].rows, splt[0].cols, splt[0].type());
        merge(splt,show);
		drawOpenCurve(show, curve, Scalar(0,0,255), 2);
#endif
		float sum_edge = 0.0;
		float sum_arc_len = 0.0;
		float arc_len_step = 5.0;
//        vector<Point2f> smooth;
//        vector<pair<char,int> > extrema = CurvatureExtrema(curve, smooth, 2.0);
//		for (int j=0; j<extrema.size()-1; j++) {
//            int i = extrema[j].second;
        for(int i=0;i<curve.size();i+=curve_step) {
        
			if (curve[i].y < dx.rows && curve[i].y >= 0 && curve[i].x < dx.cols && curve[i].x >= 0) {
//				sum_edge += dx.at<float>(curve[i].y,curve[i].x);
                
				Vec2f dir;
                if(i == curve.size()-1)
                    dir = Vec2f(curve[i].y-curve[i-1].y,-(curve[i].x-curve[i-1].x));
                else
                    dir = Vec2f(curve[i].y-curve[i+1].y,-(curve[i].x-curve[i+1].x));
				float dir_l = norm(dir);
                
                /*
				sum_arc_len += dir_l;
				if (sum_arc_len < arc_len_step) {
					continue;
				} else {
					sum_arc_len -= arc_len_step;
				}
                */

				
				dir /= dir_l;//normalize
#ifdef _DEBUG_
				line(show, curve[i], Point2f(curve[i].x+dir[0]*5,curve[i].y+dir[1]*5), Scalar(0,255), 1);
#endif
				vector<pair<float,Point2f> > edges; 
				for (int step = low; step < high; step++) {
					Point2f ev_pt(curve[i].x+dir[0]*step,curve[i].y+dir[1]*step);
                    if (ev_pt.y < dx.rows && ev_pt.y >= 0 && ev_pt.x < dx.cols && ev_pt.x >= 0) {
                        Vec3f e = dx.at<Vec3f>(cvFloor(ev_pt.y),cvFloor(ev_pt.x));
                        if (e[1] != 0 || e[2] != 0) {
                            Vec2f edge(-e[1],-e[2]);
                            float dp = edge.dot(dir) + 10 * (1.0 - (fabsf(step*2)/((float)steps)));
                            edges.push_back(make_pair(dp, ev_pt));
#ifdef _DEBUG_
                            line(show, ev_pt, Point2f(ev_pt.x+edge[0],ev_pt.y+edge[1]), Scalar(255*dp), 1);
#endif
                        }
                    }
				}
                bool penalize = false;
				if (edges.size() > 0) {
//					std::sort(edges.begin(), edges.end(), sortbyfirst<float,Point2f>);
//					if (edges.back().first > 1 && edges.back().first < 1000) {
                    vector<pair<float,Point2f> >::iterator maxedge_it = std::max_element(edges.begin(), edges.end(), sortbyfirst<float,Point2f>);
                    if ((*maxedge_it).first > 1 && (*maxedge_it).first < 1000) {
#ifdef _DEBUG_
//						circle(show, (*maxedge_it).second, 1.0+(*maxedge_it).first, Scalar(255,255,0), 1);
#endif
						edge_history.push_back((*maxedge_it).first);
//						sum_edge += (*maxedge_it).first;
						{
							p.dx_dy.push_back(make_pair(i, (*maxedge_it).second));
						}
					}
                    else {
                        penalize = true;
                    }
				} //else: penalize?
				else {
                    penalize = true;
				}
                if(penalize) {
#ifdef _DEBUG_
                    circle(show, curve[i], 5, Scalar(0,0,255), CV_FILLED);
#endif
                    edge_history.push_back(-5);
//                    sum_edge -= 15.0;
                }
			}
		}
//        //take the inliers (based on median)
//        vector<float> _eh = edge_history;
//        sort(_eh.begin(),_eh.end());
//        ShowMathMLCurve(_eh, "edge history");
//        float median = _eh[_eh.size()/2];
//        sum_edge = std::accumulate(upper_bound(_eh.begin(),_eh.end(),median),_eh.end(),0.0f);
		Scalar eh_mean,eh_stdev;
		meanStdDev(edge_history, eh_mean, eh_stdev);
        sum_edge = eh_mean[0] - eh_stdev[0];
#ifdef _DEBUG_
//		sum_edge /= edge_history.size();
		stringstream ss; ss << "sum edge " << sum_edge;

		cout << "eh_mean " << eh_mean << ", eh_stdev " << eh_stdev << "\n";
		
		//count skin pixels
        if(!skin_mask.empty()) {
            Mat msk(skin_mask.size(),CV_8UC1); fillCurve(msk, curve, Scalar::all(255));
            sum_edge += (float)(countNonZero(skin_mask & msk))/500.0;
        }

		ss << " with pixel-color " << sum_edge;
		putText(show, ss.str(), cv::Point(10,15), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(0,0,255), 1);
		imshow_<640,480>("fitness",show);
		waitKey(0);
#endif
		
		return MAX(0,sum_edge);
	}
private:
    
	/* from "Real time Hand Tracking and Gesture Recognition Using Smart Snakes", Heap & Samaria 1995 
	   section 2.2
	 */
	void localOptimizeParticles(const Mat& dx, const Mat& skin_mask) {
		vector<float> normalized_weights_ = normalized_weights;
		sort(normalized_weights_.begin(),normalized_weights_.end());
		float N_precentile = normalized_weights_[normalized_weights_.size()*0.5];

#undef _DEBUG_
#ifdef _DEBUG_
		vector<Mat> splt; split(dx, splt);
		Mat gray8; splt[0].convertTo(gray8, CV_8UC1, 255);
		
		Mat show(480,640,CV_8UC3);
#endif
#pragma omp parallel for
		for (int i=0; i<num_particles; i++) {
			//pick only strong particles
			if (normalized_weights[i] > N_precentile) 
			{
				if (particles[i].dx_dy.size() <= 0) {
					continue;
				}
				int num_iters = 500*normalized_weights[i];
                vector<Point2f> curve;
				for (int iter=0; iter < num_iters; iter++) {
#ifdef _DEBUG_
					cvtColor(gray8, show, CV_GRAY2BGR); 
#endif
					
                    if(curve.size() <= 0)
                        curve = getShapeForParticle(particles[i]);
                    
					if (particles[i].dx_dy.size() <= 0) {
                        num_iters = 0;
						break;
					}
#ifdef _DEBUG_
					drawOpenCurve(show, curve, Scalar(255), 2);
#endif
					
					vector<Point2f> original,movement;
					for (int j=0; j<particles[i].dx_dy.size(); j++) {
						original.push_back(curve[particles[i].dx_dy[j].first]);
						movement.push_back(particles[i].dx_dy[j].second);
#ifdef _DEBUG_
						circle(show, particles[i].dx_dy[j].second, 2, Scalar(0,0,255), 1);
						line(show, original.back(), movement.back(), Scalar(0,0,255), 1);
#endif
					}
					Point2f diff; float angle,scale;
					assert(original.size() > 0); assert(movement.size() > 0);
					Mat Trans = Find2DRigidTransform(original, movement, &diff, &angle, &scale);
					//transform(curve,curve,Trans);
                    
					particles[i].centroid += diff;
					particles[i].scale *= scale;
					particles[i].orientation -= (angle/CV_PI) * 180.0;
                    
					unnormalized_weights[i] = calculateFitness(particles[i], dx, skin_mask, curve, 1, norm(diff)*.75);

#ifdef _DEBUG_
					stringstream ss; ss << iter << "/" << num_iters << "("<<unnormalized_weights[i]<<")";
					putText(show, ss.str(), cv::Point(10,15), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(255), 1);
                    stringstream ss1; ss1 << "diff " << diff << "("<<norm(diff)<<")";
					putText(show, ss1.str(), cv::Point(10,25), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(255), 1);
                    stringstream ss2; ss2 << "angle " << angle;
					putText(show, ss2.str(), cv::Point(10,35), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(255), 1);
                    stringstream ss3; ss3 << "scale " << scale;
					putText(show, ss3.str(), cv::Point(10,45), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(255), 1);

					imshow_<640,480>("optimize", show);
//                    ostringstream ss4; ss4 << "local_optimize_" << i << "_" << setfill('0') << setw(5) << iter << ".png";
//                    imwrite(ss1.str(), show);
					waitKey(0);
#endif

                    if (!(norm(diff) > 2 || scale > 1.02 || scale < 0.98 || angle > 0.05 || angle < -0.05)) {
                        break;
                    }

				}
 
                if(num_iters > 0) {
//                    for(int i=0;i<2;i++) {
#undef _DEBUG_
#ifdef _DEBUG_
                    vector<Mat> splt; split(dx, splt);
                    Mat gray8; splt[0].convertTo(gray8, CV_8UC1, 255);
                    Mat grayRGB; cvtColor(gray8, grayRGB, CV_GRAY2BGR);
                    Mat show(gray8.size(),CV_8UC3);
//                    show.setTo(0);
                    double upscale = 3.33;
                    resize(grayRGB, show, Size(),upscale,upscale);
#endif
                    
                    //                Particle pc();
                    if(curve.size() <= 0)
                        curve = getShapeForParticle(particles[i]);
                    
#ifdef _DEBUG_
                    {
                        vector<Point2f> tmp(curve);
                        Mat(tmp) *= upscale;
                        drawOpenCurve(show, tmp, Scalar(255,0,0), 2);
                    }
#endif
                    //                unnormalized_weights[i] = calculateFitness(particles[i], dx, skin_mask, 1);
                    
                    for (int j=0; j<particles[i].dx_dy.size(); j++) {
                        curve[particles[i].dx_dy[j].first] = particles[i].dx_dy[j].second;
                    }
                    //                SimpleSmoothCurve(curve, curve, 3.0, true);
#ifdef _DEBUG_
                    {
                        vector<Point2f> tmp(curve);
                        Mat(tmp) *= upscale;
                        drawOpenCurve(show, tmp, Scalar(0,255,0), 2);
                    }
#endif
                    
                    Mat_<float> Tpca = Mat_<float>::eye(3,3);
                    getRotationMatrix2D(Point2f(0,0), -particles[i].orientation, 1.0/particles[i].scale).copyTo(Tpca.rowRange(0,2));
                    Tpca = Tpca * (Mat_<float>(3,3) << 1,0,-particles[i].centroid.x, 0,1,-particles[i].centroid.y, 0,0,1);
//                cout << "rot mat from angle\n"<<getRotationMatrix2D(cp, angle * 180.0 / CV_PI, 1.0) << "\n";
//                    cout << "trans\n"<<Tpca<<"\n";
                    transform(curve,curve,Tpca.rowRange(0,2));
                    
                    Mat in_shape_space = hpdm.projectToShapeSpace(curve);
                    in_shape_space.copyTo(particles[i].deformation);

                    /*
#ifdef _DEBUG_
                    {
                        Scalar coffset(250,50);
                        vector<Point2f> tmp(curve);
                        Mat(tmp) += coffset;
                        drawOpenCurve(show, tmp, Scalar(128,64,255), 1);
                        
                        tmp = hpdm.rebuildCurve(hpdm.getPatchMean(particles[i].patch_num));
                        Mat(tmp) += coffset;
                        drawOpenCurve(show, tmp, Scalar(255,64,128), 1);
                        
                        tmp = hpdm.rebuildCurve(particles[i].deformation);
                        Mat(tmp) += coffset;
                        drawOpenCurve(show, tmp, Scalar(64,255,128), 1);
                    }
#endif
                     */
                    
                    //                Mat_<double> dists;
                    //                int nearest_patch = hpdm.nearestPatchWithShapeSpaceCurve(in_shape_space,dists);
                    //                particles[i].patch_num = nearest_patch;
                    
                    //                unnormalized_weights[i] = calculateFitness(particles[i], dx, skin_mask);
                    
                    //                    Mat in_patch_space = hpdm.projectToPatchSpace(in_shape_space,pc.patch_num);
                    //                    hpdm.backprojectFromPatchSpace(in_patch_space, pc.patch_num).copyTo(particles[i].deformation);
#ifdef _DEBUG_
                    curve = getShapeForParticle(particles[i]);
                    {
                        vector<Point2f> tmp(curve);
                        Mat(tmp) *= upscale;
                        drawOpenCurve(show, tmp, Scalar(0,0,255), 2);
                    }
                    imshow("active contour", show);
                    waitKey(0);
#endif
//                    curve.clear();
                    
                    unnormalized_weights[i] = calculateFitness(particles[i], dx, skin_mask, curve, 2);
//                    }
                }
                
			}
		}
		normalizeWeights();
	}
    
    /**
     * calculate the number of effective particles. Based on http://en.wikipedia.org/wiki/Particle_filter on SIS
     * @return number of effective particles
     **/
    float numEffectiveParticles() {
        //N_eff = 1 / sum(w_i)^2
        float N_eff_denom = 0.0;
        for(int i=0;i<num_particles;i++) {
            N_eff_denom += normalized_weights[i] * normalized_weights[i];
        }
        return 1.0 / N_eff_denom;
    }
	
	void normalizeWeights() {
		// from http://en.wikipedia.org/wiki/Particle_filter on SIR:
		// w_t(i) = w_t-1(i) * p(z|x);
	
//        vector<float> last_unnorm(unnormalized_weights);
//		double maxv; minMaxLoc(unnormalized_weights, 0, &maxv);
		for (int i=0; i<unnormalized_weights.size(); i++) {
			//minimizing
//			unnormalized_weights[i] = normalized_weights[i] * (1 + maxv - unnormalized_weights[i]);
			//maximizing
			unnormalized_weights[i] = (1e-2 + normalized_weights[i]) * unnormalized_weights[i];
		}
	
		
		//sum all weights
		Scalar s = sum(unnormalized_weights);
		
		//calc normalized weights by dividing, and build cumulative distribution function
        cumulative_sum[0] = 0.0;
		for (int i=0; i<unnormalized_weights.size(); i++) {
			normalized_weights[i] = unnormalized_weights[i] / s[0];
//			if (i==0) {
//				cumulative_sum[0] = normalized_weights[0];
//			} else {
            if(i>0)
				cumulative_sum[i] = cumulative_sum[i-1] + normalized_weights[i];
//			}
		}
	}	
	
	void calculateFitnessForAllParticles(const Mat& dx, const Mat& dy) {
#pragma omp parallel for
		for (int i=0; i<num_particles; i++) {
            vector<Point2f> curve;
			unnormalized_weights[i] = calculateFitness(particles[i], dx, dy, curve);
//			cout << "particle " <<i<< " " <<unnormalized_weights[i]<<"\n";
		}
	}
	
	void add_noise(Particle& p) {
		//add noise based on associated patch (subregion)
		//major axes for subregion
//		vector<Point2f> curve = getShapeForParticle(p);
//		Mat curve_as_mat = Mat(curve).reshape(1,curve.size());
//		PCA cpca(curve_as_mat,Mat(),CV_PCA_DATA_AS_ROW);
//		cout << "pca for particle \n" <<cpca.eigenvectors<<"\n"<<cpca.eigenvalues<<"\n";
		//add position noise based on major axes
//		Mat axis1 = cpca.eigenvectors.col(0) * rng_.gaussian(sqrt(cpca.eigenvalues.at<float>(0))/8.0);
//		Mat axis2 = cpca.eigenvectors.col(1) * rng_.gaussian(sqrt(cpca.eigenvalues.at<float>(1))/8.0);
//		Mat movement = axis1+axis2;
//		p.centroid += Point2f(movement.at<float>(0),movement.at<float>(1));
        
        //TODO: all these parameters can be learned from a clean run of the tracker, or use a KF
        p.centroid += Point2f(rng_.gaussian(1.0),rng_.gaussian(1.0));
		p.orientation += rng_.gaussian(3.0);
		p.scale += rng_.gaussian(0.05);
        
        
//		float pa_sigma = 0;
//		Mat deformation_noise = hpdm.getPatchPrincipalAxis(p.patch_num,pa_sigma);
//////		rng_.fill(deformation_noise,RNG::NORMAL,0.0,0.05);
////		//TODO: deformation noise should be distributed by the patch PCA
//////		p.deformation += deformation_noise;
//		addWeighted(p.deformation,1.0,deformation_noise,rng_.gaussian(pa_sigma),0.0,p.deformation);
	}
	
	void resampleParticles() {
		//draw N particles from the new distribution function (the importance density)
		vector<Particle> new_particles(num_particles);
		for (int i=0; i<particles.size(); i++) {
			while (true) {
				float uni_rand = rng_.uniform(0.0f,1.0f);
				//binary search for position in cumulative sum
				vector<float>::iterator pos = lower_bound(cumulative_sum.begin(), cumulative_sum.end(), uni_rand);
				int ipos = distance(cumulative_sum.begin(), pos);
                
                //re-sample if scale is too big
                if(particles[ipos].scale < 0.75*initscale || particles[ipos].scale > 1.25*initscale) continue;
				
                new_particles[i].scale = particles[ipos].scale;
				new_particles[i].orientation = particles[ipos].orientation;
                
                //re-sample if out of image bounds
                if(particles[ipos].centroid.x < 0 || particles[ipos].centroid.y < 0 ||
                   particles[ipos].centroid.x > bounds.width || particles[ipos].centroid.y > bounds.height)
                    continue;
                
				new_particles[i].centroid = particles[ipos].centroid;
				new_particles[i].patch_num = particles[ipos].patch_num;
				particles[ipos].deformation.copyTo(new_particles[i].deformation);
				break;
			}
		}
		
		particles.clear(); particles.insert(particles.begin(),new_particles.begin(),new_particles.end());
        for (int i=0; i<num_particles; i++) {
            normalized_weights[i] = 1.0 / (float)num_particles;
        }
	}
	
	void propagate() {
		//according to section 3.2
		int changed = 0;
		for (int i=0; i<particles.size(); i++) {
			int old_patch_num = particles[i].patch_num;
            
			float uni_rand = rng_.uniform(0.0f,1.0f);
            
			//binary search for new subregion (patch) in Cumulative Transition Matrix
            double* cumulative_ptr = Cumulative.ptr<double>(old_patch_num);
            vector<double> row_a_in_CumulativeMat(cumulative_ptr,cumulative_ptr+Cumulative.cols);
			vector<double>::iterator pos = upper_bound(row_a_in_CumulativeMat.begin(), row_a_in_CumulativeMat.end(), uni_rand);
			int new_patch_num = distance(row_a_in_CumulativeMat.begin(), pos);

            
            //			cout << "particle "<<i<<" ("<<uni_rand<<") ["<<particles[i].patch_num<<"] -> ["<<new_patch_num<<"]\n";
			if (new_patch_num != old_patch_num) {
                if(new_patch_num >= 0 && new_patch_num < hpdm.getNumPatches()) {
                    //moved to a new patch ("leap through wormhole"): take the mean deformation of that patch
                    particles[i].patch_num = new_patch_num;
                    hpdm.getPatchMean(new_patch_num).copyTo(particles[i].deformation);
                    changed++;
                }
			}
             
			 
			add_noise(particles[i]);
			
			//TODO: particle should not go out of bounds or be too small/big or orientate too much
		}
//		cout << changed << " particles changed, " << (num_particles-changed) << " stayed\n";
	}
		
	void calculateMeanShape() {
        //Take the particles that score the highest
		vector<float> normalized_weights_ = normalized_weights;
		sort(normalized_weights_.begin(),normalized_weights_.end());
		float N_precentile = normalized_weights_[normalized_weights_.size()*0.95];
//        float N_precentile = *(max_element(normalized_weights.begin(),normalized_weights.end())) * 0.95;
        if(isnan(N_precentile)) return;

        mean_particle.scale = mean_particle.orientation = mean_particle.centroid.x = mean_particle.centroid.y = 0;
        if(mean_particle.deformation.empty())
            mean_particle.deformation.create(1,hpdm.getShapeSpaceDim());
		mean_particle.deformation.setTo(0);

		float sum_w = 0.0;
		for (int i=0; i<num_particles; i++) {
			if (normalized_weights[i] < N_precentile) {
				continue;
			}
			mean_particle.orientation += particles[i].orientation * normalized_weights[i];
			mean_particle.scale += particles[i].scale * normalized_weights[i];
			mean_particle.centroid += particles[i].centroid * normalized_weights[i];
			addWeighted(mean_particle.deformation, 1.0, particles[i].deformation, normalized_weights[i], 0.0, mean_particle.deformation);
			sum_w += normalized_weights[i];
		}
		mean_particle.orientation /= sum_w;
		mean_particle.scale /= sum_w;
		mean_particle.centroid *= 1.0/sum_w;
		mean_particle.deformation /= sum_w;

		mean_shape = getShapeForParticle(mean_particle);
	}
    
    void initializeInternal() {
        assert(hpdm.isInitialized());
        
        m_curves = hpdm.getAlignedCurves();
        
		buildTransitionMatrix();
		
		initializeParticles();
        
		this->initialized = true;
	}

	
public:
	HPDMTracker():num_particles(30),curve_length(120),m_iter(0),initscale(1.0) {
        rng_ = RNG(time(NULL));
        N_thr = num_particles*0.25;
    }
		
    vector<Point2f> getShapeForParticle(const Particle& p) const {
		vector<Point2f> shape = hpdm.rebuildCurve(p.deformation);
		//transform according to particle
		Mat_<float> Tpca = Mat_<float>::eye(3,3);
		getRotationMatrix2D(Point2f(0,0), p.orientation, p.scale).copyTo(Tpca.rowRange(0,2));
		Tpca = (Mat_<float>(3,3) << 1,0,p.centroid.x, 0,1,p.centroid.y, 0,0,1) * Tpca;
		//	cout << "rot mat from angle\n"<<getRotationMatrix2D(cp, angle * 180.0 / CV_PI, 1.0) << "\n";
		//		cout << "trans\n"<<Tpca<<"\n";
		transform(shape,shape,Tpca.rowRange(0,2));
		return shape;
	}

    void saveToFile(const char* filename) {
        hpdm.saveToFile(filename);
    }
        
	void initializeFromDirectory(const char* dir_name) {
		vector<string> curves_files = open_dir(dir_name, "txt");
				
		Mat show(480,640,CV_8UC3);
		
		cout << "reading from file...\n";
		vector<vector<Point2f> > curves;
		Mat curves_by_row(curves_files.size(),curve_length*2,CV_32FC1);
		for (int i=0; i<curves_files.size(); i++) {            
			vector<Point2f> c_ = loadCurveFromFile<float>(curves_files[i]);
			vector<Point2f> c; c.insert(c.begin(),c_.begin()+120,c_.end()-100);
			ResampleCurve(c, c, curve_length); //normalize to N points
			curves.push_back(c);
			Mat(Mat(c).t()).reshape(1, 1).copyTo(curves_by_row.row(i));
			
			drawOpenCurve(show, c, Scalar(0,0,255), 1);
		}
		vector<Point2f> avg_curve_v(curve_length);
		Mat curves_mean;
		reduce(curves_by_row, curves_mean, 0, CV_REDUCE_AVG);
		Mat(curves_mean.t()).reshape(2, curve_length).copyTo(Mat(avg_curve_v));	
		drawOpenCurve(show, avg_curve_v, Scalar(255), 2);
		
//		imshow("original DB", show);
//		waitKey(1);
		
		initialize(curves);
	}		
	
    void initializeFromFile(const char* filename, bool realign) {
        hpdm.initializeFromFile(filename);
        if(!hpdm.isInitialized()) {
            cerr << "cannot initialize HPDM" << endl;
            return;
        }
        if(realign)
            hpdm.realign();
        curve_length = hpdm.getCurveLength();
        initializeInternal();
    }

	void initialize(const vector<vector<Point2f> >& sequential_curves) {
		hpdm.initialize(sequential_curves);
		destroyAllWindows();
        
        initializeInternal();
    }
    
	void setInitialGuess(cv::Point center, float scale, float orientation) {
        initscale = scale;
		for (int i=0; i<num_particles; i++) {
			//randomize based on guess
			particles[i].orientation = orientation; // + rng_.gaussian(5);
			particles[i].centroid.x = center.x + rng_.gaussian(.5);
			particles[i].centroid.y = center.y + rng_.gaussian(.5);
			particles[i].scale = scale + rng_.gaussian(0.1);
			particles[i].patch_num = rng_.uniform(0, hpdm.getNumPatches()-1);
			hpdm.getPatchMean(particles[i].patch_num).copyTo(particles[i].deformation);
            normalized_weights[i] = 1.0 / (float)num_particles;
		}
	}
    void setInitialGuess(cv::Point center, float scale, float orientation, int patch_num) {
        setInitialGuess(center, scale, orientation);
		for (int i=0; i<num_particles; i++) {
			particles[i].patch_num = patch_num;
			hpdm.getPatchMean(particles[i].patch_num).copyTo(particles[i].deformation);
		}
    }
    
    void visualizeParticles(Mat& show) const {
        cv::Point high_loc;
        double minVal,maxVal;
        //         CV_PROFILE_MSG("show",
        minMaxLoc(normalized_weights, &minVal, &maxVal,0,&high_loc);
        for (int i=0; i<num_particles; i++) {
            vector<Point2f> s = getShapeForParticle(particles[i]);
            if (i==high_loc.x) {
                drawOpenCurve(show, s, Scalar(0,255,255), 2);
            } else {
                drawOpenCurve(show, s, Scalar(0,0,255*(normalized_weights[i]-minVal)/(maxVal-minVal)),1);
            }
        }
        //		calculateMeanShape();
        drawOpenCurve(show, mean_shape, Scalar(255,0,255), 2);
        //         imshow("particles", show);
        //                       stringstream filename; filename << "pf_" << setfill('0') << setw(4) << m_iter << ".jpg";
        //                       imwrite(filename.str(), show);
        //         )
    }
	
    Mat gray;
    Mat gray32f;
    Mat dx;
	void update(const Mat& rgb_img, const Mat& skin_mask) {
        bounds = rgb_img.size();
        
		CV_PROFILE_MSG("prepare",
        if(rgb_img.channels() == 3)
            cvtColor(rgb_img, gray, CV_BGR2GRAY);
        else
            gray = rgb_img;
//		equalizeHist(gray, gray);		
        gray.convertTo(gray32f, CV_32FC1, 1.0/255.0);
//		Mat gray32f; skin_mask.convertTo(gray32f, CV_32FC1, 1.0/255.0);
		vector<Mat> dv(3);
        gray32f.copyTo(dv[0]);
        Sobel(gray32f, dv[1], -1, 1, 0, CV_SCHARR);
        Sobel(gray32f, dv[2], -1, 0, 1, CV_SCHARR);
//		double minVal,maxVal;
//		minMaxLoc(dx, &minVal, &maxVal);// cout << "minval " << minVal << " maxval " << maxVal << "\n";
//		dx = abs(dx)/maxVal;
//		minMaxLoc(dy, &minVal, &maxVal);
//		dy = abs(dy)/maxVal;
//		dx = (dx + dy) / 2.0;
//		imshow("gray32f", gray32f);
		merge(dv,dx);
		
//		imshow("edge",dx);
//		imshow("skin",skin_mask);
//		waitKey(1);
				   )
		
		CV_PROFILE(calculateFitnessForAllParticles(dx,skin_mask);)
		CV_PROFILE(normalizeWeights();)
		
		CV_PROFILE(localOptimizeParticles(dx,skin_mask);)
//		CV_PROFILE(normalizeWeights();)

		CV_PROFILE(calculateMeanShape();)

        
         
//		Mat curve_as_mat = Mat(mean_shape).reshape(1,mean_shape.size());
//		PCA cpca(curve_as_mat,Mat(),CV_PCA_DATA_AS_ROW);
//		Point c(cpca.mean.at<float>(0),cpca.mean.at<float>(1));
//		Mat axis1 = cpca.eigenvectors.col(0) * sqrt(cpca.eigenvalues.at<float>(0));
//		Mat axis2 = cpca.eigenvectors.col(1) * sqrt(cpca.eigenvalues.at<float>(1));
//		line(show, c, c+Point(axis1.at<float>(0),axis1.at<float>(1)), Scalar(0,255), 2);
//		line(show, c, c+Point(axis2.at<float>(0),axis2.at<float>(1)), Scalar(255), 2);
		
        N_eff = numEffectiveParticles();
        if(isnan(N_eff)) N_eff = 0.0f;
        if(N_eff < N_thr) {
            cout << "N_eff "<<N_eff<<" < "<<N_thr<<" (N_thr) :";
//            int N_eff_i = (int)N_eff;
//            ShowMathGLDataAsHist(normalized_weights,&(N_eff_i)/*,"before resample"*/);
            if(N_eff <= 0) {
                //if N_eff is 0 it means no particle is effective, better default to last known position
                cout << " re-init from last known particle\n";
                setInitialGuess(mean_particle.centroid, mean_particle.scale, mean_particle.orientation);
            } else {
                cout << " resample\n";
                CV_PROFILE(resampleParticles();)
            }
//            ShowMathGLDataAsHist(normalized_weights,NULL,"after resmaple");
//            waitKey();
        }
//        else {
//            cout << "N_eff "<<N_eff<<endl;
//        }
        
		CV_PROFILE(propagate();)
        
        m_iter++;
	}
	
	vector<Point2f> getMeanShape() const { return mean_shape; }
    int getMeanShapePatch() const {
        //return the patch of the highest ranking particle
        return particles[distance(normalized_weights.begin(), max_element(normalized_weights.begin(), normalized_weights.end()))].patch_num;
    }
    const HPDM<float>& getHPDM() const { return hpdm; }
    const float getNumEffectiveParticles() const { return N_eff; }
};