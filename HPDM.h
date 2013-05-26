/*
 *  HPDM.h
 *  CurveMatching
 *
 *  Created by roy_shilkrot on 2/15/13.
 *  Copyright 2013 MIT. All rights reserved.
 *
 * from "Improving Specificity in PDMs using Hierarchical Approach", Heap and Hogg, 1997
 * and "Wormholes in Shape Space: Tracking through Discontinuous Changes in Shape", Heap & Hogg, 1998
 */
#pragma once

#include "AbstractAlgorithm.h"

template<typename T>
class HPDM : public AbstractAlgorithm {
	
	vector<vector<Point_<T> > > original_curves;
	vector<vector<Point_<T> > > aligned_curves;
	Mat aligned_curves_by_row;
	int curve_length;
	int shape_space_dim;
    int num_pc_for_subregions;
	vector<Point2f> curves_mean_v;
	
	PCA curves_pca;
	vector<PCA> subregions_pca;
	Mat subregion_means;
	vector<T> subregions_lambda;
    
    void addRemovePointsToMatchCurve(const vector<Point2f>& mean_curve_orig, vector<Point2f>& curve) {
        if(norm(curve[0] - mean_curve_orig[0]) > 2) {
            int closest = 0; double minD = norm(curve[0] - mean_curve_orig[0]), newD = -1;
            while ((newD = norm(curve[0] - mean_curve_orig[closest])) <= minD) {
                closest++;
                minD = newD;
            }
            //                cout << "start closest add " << closest << " mind " << minD << "\n";
            //make up for the space by adding a piece from the last curve
            if (closest > 1) {
                //                    circle(show, original_curves[i][0], 2, Scalar(0,255),CV_FILLED);
                Point2f off = curve[0] - mean_curve_orig[closest];
                for (int j=closest-1; j>=0; j--) {
                    curve.insert(curve.begin(),mean_curve_orig[j] + off);
                }
                //                    circle(show, original_curves[i][0], 2, Scalar(255,255),CV_FILLED);
            } else { //nothing added, maybe remove?
                closest = 0;
                minD = norm(curve[0] - mean_curve_orig[0]);
                while ((newD = norm(curve[closest] - mean_curve_orig[0])) <= minD) {
                    closest++;
                    minD = newD;
                }
                //                    cout << "start closest erase " << closest << " mind " << minD << "\n";
                curve.erase(curve.begin(),curve.begin()+(closest-1));
            }
        }
        
        //check if the curves end roughly at the same point
        if(norm(curve.back() - mean_curve_orig.back()) > 2) {
            int closest = mean_curve_orig.size()-1;
            double minD = norm(curve.back() - mean_curve_orig.back()), newD = -1;
            while ((newD = norm(curve.back() - mean_curve_orig[closest])) <= minD) {
                closest--;
                minD = newD;
            }
//                            cout << "end closest add " << closest << " mind " << minD << "\n";
            if (closest < mean_curve_orig.size()-2) {
                Point2f off = curve.back() - mean_curve_orig[closest];
                for (int j=closest+1; j<mean_curve_orig.size(); j++) {
                    curve.push_back(mean_curve_orig[j] + off);
                }
            } else { //nothing added, maybe remove?
                closest = curve.size()-1;
                minD = norm(curve.back() - mean_curve_orig.back());
                while ((newD = norm(curve[closest] - mean_curve_orig.back())) <= minD) {
                    closest--;
                    minD = newD;
                }
//                                   cout << "end closest erase " << closest << " mind " << minD << "\n";
                curve.erase(curve.begin()+closest+1,curve.end());
            }
            
        }

    }
    
    struct HPDM_onMouse_data {
        void* curve;
        bool mouseAdvance;
        int idx;
        int tot;
    };

    static void HPDM_onMouse(int e, int x, int y, int flags, void* data_) {
        static Mat show;
        
        if(show.empty())
            show.create(500,500,CV_8UC3);

        HPDM_onMouse_data* data = static_cast<HPDM_onMouse_data*>(data_);
        assert(data);
        
        vector<Point2f>* curve_v = static_cast<vector<Point2f>* >(data->curve);
        assert(curve_v);
        
        show.setTo(0);
        vector<Point2f> tmp; ConvertCurve(*curve_v,tmp);
        Mat(tmp) += Scalar(150,150);
        drawOpenCurve(show, tmp, Scalar(0,0,255), 1);
        
        int minLoc = closestPointOnCurveToPoint(tmp, cv::Point(x,y), 15.0f);
        if(minLoc > -1) {
            line(show, cv::Point(x,y), tmp[minLoc], Scalar(0,255));
            circle(show, tmp[minLoc], 4, Scalar(255), CV_FILLED);
        }
        
        if( e == CV_EVENT_LBUTTONDOWN ) {
            //press
        } else if( e == CV_EVENT_LBUTTONUP ) {
            //mouse release
            if(minLoc > -1) {
                Mat(tmp) -= Scalar(150,150);
                if(minLoc > (tmp.size() - minLoc)) {
                    //closer to the end
                    tmp.erase(tmp.begin() + minLoc,tmp.end());
                } else {
                    //closet to begining
                    tmp.erase(tmp.begin(),tmp.begin() + minLoc);
                }
                ResampleCurve(tmp, tmp, curve_v->size());
                curve_v->clear(); curve_v->insert(curve_v->begin(),tmp.begin(),tmp.end());
                data->mouseAdvance = true;
            }
        } else if (e == CV_EVENT_MOUSEMOVE && flags == CV_EVENT_FLAG_LBUTTON) {
            //drag
        }

        stringstream ss; ss << data->idx << "/" << data->tot;
        putText(show,ss.str(),cv::Point(10,15),CV_FONT_HERSHEY_PLAIN,1.0,Scalar::all(255));
        imshow("aligned", show);
    }

	
	void align() {
		/*********** HH 1997 Section 2: Aligning ***********/
		
		cout << "aligning...\n";
		//Align mean curve
        int first_curve = 40;
		vector<Point2f> mean_curve_orig = original_curves[first_curve];
        vector<Point2f> mean_curve = mean_curve_orig;
		//	cout << "Mat(mean_curve).reshape(1,curve_length)\n " << Mat(mean_curve).reshape(1,curve_length) << "\n";
		PCA mean_curve_pca(Mat(mean_curve).reshape(1,curve_length),Mat(),CV_PCA_DATA_AS_ROW);
		
        cout << "mean curve pca ev\n"<<mean_curve_pca.eigenvectors<<"\n";
		
        float angle = atan2(mean_curve_pca.eigenvectors.at<float>(0,1), mean_curve_pca.eigenvectors.at<float>(0,0));
		float scale = sqrt(mean_curve_pca.eigenvalues.at<float>(0,0));
		Point2f cp = Point2f(mean_curve_pca.mean.at<float>(0),mean_curve_pca.mean.at<float>(1));
        
		cout << "cp " << cp << " angle " << angle << " scale " << scale << "\n";
		
        Mat_<float> Tpca(3,3);
		getRotationMatrix2D(Point2f(0,0), angle * 180.0 / CV_PI, 1.0).copyTo(Tpca.rowRange(0,2));
		Tpca *= (Mat_<float>(3,3) << 1,0,-cp.x, 0,1,-cp.y, 0,0,1);
		//	cout << "rot mat from angle\n"<<getRotationMatrix2D(cp, angle * 180.0 / CV_PI, 1.0) << "\n";
		//	Mat(mean_curve_pca.eigenvectors.t()).copyTo(Tpca.colRange(0, 2));
		
        cout << "trans\n"<<Tpca<<"\n";
		
        transform(mean_curve,mean_curve,Tpca.rowRange(0,2));
		
		Mat show(500,500,CV_8UC3,Scalar::all(0));
		//	drawOpenCurve(show, mean_curve, Scalar(255), 2);
		
		aligned_curves_by_row.create(original_curves.size(),curve_length*2, CV_32FC1);
		for (int i=0; i<original_curves.size(); i++) {
            if (i==first_curve) {
                continue;
            }
//            show.setTo(0);
//			drawOpenCurve(show, mean_curve_orig, Scalar(255,255,255), 1);
//            circle(show,mean_curve_orig[100],5,Scalar(128,255,0),1);
//			drawOpenCurve(show, original_curves[i], Scalar(255,0,255), 1);
//            circle(show,original_curves[i][100],5,Scalar(255,0,255),1);
            
            //check if the curves start roughly at the same point
            Scalar mean_a = mean(original_curves[i]),mean_b = mean(mean_curve_orig);
            Mat_<float> move_to_mean = (Mat_<float>(2,3) << 1,0,-mean_a[0]+mean_b[0],  0,1,-mean_a[1]+mean_b[1]);
            transform(original_curves[i], original_curves[i], move_to_mean);
            
            addRemovePointsToMatchCurve(mean_curve_orig,original_curves[i]);
            
            ResampleCurve(original_curves[i], original_curves[i], mean_curve.size());
// 			drawOpenCurve(show, original_curves[i], Scalar(0,0,255), 1);
//            circle(show,original_curves[i][100],5,Scalar(0,255,255),CV_FILLED);
           
//			Mat Trans = Find2DRigidTransform(original_curves[i], mean_curve);
//            //            cout << "trans\n"<<Trans<<"\n";
			vector<Point2f> transformed(mean_curve.size());
//			transform(original_curves[i],transformed,Trans);
//            aligned_curves.push_back(transformed);


            Mat Trans = Find2DRigidTransform(original_curves[i], mean_curve_orig);
			transform(original_curves[i],transformed,Trans);
            
//			drawOpenCurve(show, transformed, Scalar(0,0,255), 1);
//            circle(show,transformed[100],5,Scalar(0,255,255),CV_FILLED);
            
            addRemovePointsToMatchCurve(mean_curve_orig,transformed);
            ResampleCurve(transformed, transformed, mean_curve.size());
            
//            drawOpenCurve(show, transformed, Scalar(0,0,255), 1);
//            circle(show,transformed[100],5,Scalar(0,255,255),CV_FILLED);
            
            Trans = Find2DRigidTransform(transformed, mean_curve);
			transform(transformed,transformed,Trans);

            show.setTo(0);
            {
                vector<Point2f> transformed_tmp; ConvertCurve(transformed,transformed_tmp);
                Mat(transformed_tmp) += Scalar(150,150);
                drawOpenCurve(show, transformed_tmp, Scalar(0,0,255), 1);
                stringstream ss; ss << i << "/" << original_curves.size();
                putText(show,ss.str(),cv::Point(10,15),CV_FONT_HERSHEY_PLAIN,1.0,Scalar::all(255));
            }
            
            namedWindow("aligned");
            HPDM_onMouse_data data = {&transformed,false,i,original_curves.size()};
            setMouseCallback("aligned", HPDM::HPDM_onMouse, &data);
            imshow("aligned", show);
            while(!data.mouseAdvance)
                if(waitKey(1)==' ') break;

            aligned_curves.push_back(transformed);
            
            Mat transformed_mat = Mat(transformed).t();
            Mat transformed_mat_one_ch = transformed_mat.reshape(1, 1);
			transformed_mat_one_ch.copyTo(aligned_curves_by_row.row(i));
        }
		
		Mat curves_mean;
		curves_mean_v.clear(); curves_mean_v.resize(curve_length);
		reduce(aligned_curves_by_row, curves_mean, 0, CV_REDUCE_AVG);
		Mat(curves_mean.t()).reshape(2,curve_length).copyTo(Mat(curves_mean_v));
		
//        {
//        Scalar mean_curve_orig_mean = mean(mean_curve_orig);
//        Mat_<float> move_to_mean = (Mat_<float>(2,3) << 1,0,mean_curve_orig_mean[0],  0,1,mean_curve_orig_mean[1]);
//        transform(curves_mean_v, curves_mean_v, move_to_mean);
//        }

//		drawOpenCurve(show, curves_mean_v, Scalar(255), 2);
//		imshow("aligned", show);
//		waitKey(1);
		
//        show.setTo(0);
//        for (int i=0; i<aligned_curves.size(); i++) {
//            drawOpenCurve(show, aligned_curves[i], Scalar(255), 1, cv::Point(150,150));
//        }
//        drawOpenCurve(show, curves_mean_v, Scalar(0,0,255), 1, cv::Point(150,150));
//        imshow("aligned", show);
//        waitKey();
	}
    

	
	void buildHPDM() {
        curves_pca = PCA();
        subregions_pca.clear();
        subregion_means = Mat();
        subregions_lambda.clear();

        
		/*********** HH 1997 Section 3: HPDM ***********/
		
		float O = 1.5;
		int E = original_curves.size();
        int divisor = 30;
		int k = E / divisor;
        int NNs = divisor * O;

		cout << "global PCA...\n";
		CV_PROFILE(curves_pca(aligned_curves_by_row,Mat(),CV_PCA_DATA_AS_ROW,shape_space_dim);)
		
		Mat compressed_curves(aligned_curves_by_row.rows,shape_space_dim,CV_32FC1);
		for (int i=0; i<aligned_curves_by_row.rows; i++) {
			curves_pca.project(aligned_curves_by_row.row(i),compressed_curves.row(i));
		}
		
		Mat labels,centers;
		cout << "kmeans...("<<k<<" centers)\n";
		kmeans(compressed_curves, 
			   k, 
			   labels, 
			   TermCriteria(TermCriteria::COUNT + TermCriteria::EPS,1000,0.0001), 
			   50, //try 
			   KMEANS_RANDOM_CENTERS,
			   centers);
		
        
        int squaresize = 200;
        int sq_k = cvCeil(sqrt(centers.rows));
		Mat show(squaresize*sq_k,squaresize*sq_k,CV_8UC3,Scalar::all(0));
		
        int countcenters = 0;
		Mat good_centers_as_rows(centers.rows,centers.cols,CV_32FC1);
		for (int i=0; i<centers.rows; i++) {
			Mat rebuilt;
			curves_pca.backProject(centers.row(i),rebuilt);
            
            Scalar mn,sdv;
            meanStdDev(Mat(rebuilt).reshape(2,curve_length), mn, sdv);
            cout << "center " <<  i<< " mean " <<mn << " sdv "<<sdv<<" nsdv " << norm(sdv) <<"\n";
            if(norm(sdv) < 30) { //TODO: this constant is dependent of the data...
                cerr << "this curve is problematic\n";
                continue;
            }
            
			centers.row(i).copyTo(good_centers_as_rows.row(countcenters++));
            
			vector<Point2f> rebuilt_curve(curve_length);
			Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
			Mat(rebuilt_curve) += Scalar(squaresize/2,squaresize/2);
            Mat showroi = show(Range(squaresize*(i/sq_k),squaresize*(i/sq_k + 1)),Range(squaresize*(i%sq_k),squaresize*(i%sq_k + 1)));
			drawOpenCurve(showroi, rebuilt_curve, Scalar(0,0,255), 1);
		}
//		imshow("centers",show);
//		waitKey(0);
		k = countcenters;
        centers.create(countcenters,good_centers_as_rows.cols,centers.type());
        good_centers_as_rows.rowRange(0, countcenters).copyTo(centers);
				
		//Find nearest neighbors in the dataset for each center
		vector<vector<DMatch> > matches;
		BFMatcher bfm;
		bfm.knnMatch(centers, compressed_curves, matches, NNs);
		
		vector<Mat> centers_NNs(k);
		cout << "Find cluster means with NN...\n";
		//Find cluster means
		subregion_means.create(centers.rows,shape_space_dim,CV_32FC1);
        show.setTo(0);
		for (int i=0; i<matches.size(); i++) {

			Mat rebuilt;

			centers_NNs[i].create(NNs,shape_space_dim,CV_32FC1);
			for (int j=0; j<matches[i].size(); j++) {
				compressed_curves.row(matches[i][j].trainIdx).copyTo(centers_NNs[i].row(j));
/*
				curves_pca.backProject(compressed_curves.row(matches[i][j].trainIdx),rebuilt);
				rebuilt.copyTo(rebuilt_centers_as_rows.row(i));
				vector<Point2f> rebuilt_curve(curve_length);
				Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
				Mat(rebuilt_curve) += Scalar(150,150);
				drawOpenCurve(show, rebuilt_curve, Scalar(0,255,255), 1);
*/				
			}
			reduce(centers_NNs[i], subregion_means.row(i), 0, CV_REDUCE_AVG);

            curves_pca.backProject(subregion_means.row(i),rebuilt);

            vector<Point2f> rebuilt_curve(curve_length);
            Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
			Mat(rebuilt_curve) += Scalar(squaresize/2,squaresize/2);
            Mat showroi = show(Range(squaresize*(i/sq_k),squaresize*(i/sq_k + 1)),Range(squaresize*(i%sq_k),squaresize*(i%sq_k + 1)));
			drawOpenCurve(showroi, rebuilt_curve, Scalar(0,0,255), 1);

            /*
			curves_pca.backProject(subregion_means.row(i),rebuilt);
			rebuilt.copyTo(rebuilt_centers_as_rows.row(i));
			vector<Point2f> rebuilt_curve(curve_length);
			Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
			Mat(rebuilt_curve) += Scalar(150,150);
			drawOpenCurve(show, rebuilt_curve, Scalar(0,0,255), 1);
			imshow("subregions means",show);
			waitKey(0);
             */
		}
//        imshow("subregions means",show);
//        waitKey(0);

		
		cout << "subregions pca...\n";
		//PCA subregions
		subregions_pca.resize(k);
		subregions_lambda.resize(k);
//		subregion_means.create(k,shape_space_dim,CV_32FC1);
		vector<int> to_remove;
		for (int i=0; i<k; i++) {
			
			//get samples that belong to this region
			Mat subregion_samples_m(countNonZero(labels == i),shape_space_dim,CV_32FC1);
			if (subregion_samples_m.rows < 5 /*shape_space_dim*/) {
				to_remove.push_back(i);
				cout << "remove subregion " << i << "("<<subregion_samples_m.rows<<")\n";
				continue; //subregion too small!
			}
			int count = 0;
			for (int j=0; j<labels.rows; j++) {
				if (labels.at<int>(j) == i) {
					compressed_curves.row(j).copyTo(subregion_samples_m.row(count++));
				}
			}
			 assert(count == subregion_samples_m.rows);
			 
			/*
			if (matches[i].size() < 10) {
				to_remove.push_back(i);
				cout << "remove subregion " << i << "\n";
				continue; //subregion too small!
			}
			 //use the NN-mean for PCA
			 subregions_pca[i](centers_NNs[i],subregion_means.row(i),CV_PCA_DATA_AS_ROW);
			 */

			//use the NN-mean for PCA
			subregions_pca[i](subregion_samples_m,subregion_means.row(i),CV_PCA_DATA_AS_ROW,num_pc_for_subregions);
			float sum_ev = sum(subregions_pca[i].eigenvalues)[0];
			subregions_lambda[i] = sqrt(sum_ev);
//			subregions_pca[i].mean.copyTo(subregion_means.row(i));
		}
		if (to_remove.size() > 0) {
			//remove "empty" subregions
			Mat new_subregions_means(subregion_means.rows-to_remove.size(),subregion_means.cols,subregion_means.type());
			int count = subregion_means.rows-to_remove.size()-1;
			for (int i=k-1; i >= 0; i--) {
				if (find(to_remove.begin(), to_remove.end(), i)==to_remove.end()) {
					subregion_means.row(i).copyTo(new_subregions_means.row(count--));
				} else {
					subregions_pca.erase(subregions_pca.begin() + i);
					subregions_lambda.erase(subregions_lambda.begin() + i);
				}
			}
			subregion_means = new_subregions_means;
		}
        
	}
    
    void clearDataStruct() {
        aligned_curves_by_row = Mat();
        aligned_curves.clear();
        curve_length = shape_space_dim = num_pc_for_subregions = -1;
        curves_mean_v.clear();
        curves_pca = PCA();
        subregions_pca.clear();
        subregion_means = Mat();
        subregions_lambda.clear();
    }
    
public:
	HPDM():shape_space_dim(20),num_pc_for_subregions(2) {}
	
    void realign() {
        cout << "re-aligning...\n";
		//Align mean curve
        int first_curve = 1;
		vector<Point2f> mean_curve_orig = aligned_curves[first_curve];
        vector<Point2f> mean_curve = mean_curve_orig;
		//	cout << "Mat(mean_curve).reshape(1,curve_length)\n " << Mat(mean_curve).reshape(1,curve_length) << "\n";
		PCA mean_curve_pca(Mat(mean_curve).reshape(1,curve_length),Mat(),CV_PCA_DATA_AS_ROW);
		
        cout << "mean curve pca ev\n"<<mean_curve_pca.eigenvectors<<"\n";
		
        float angle = atan2(mean_curve_pca.eigenvectors.at<float>(0,1), mean_curve_pca.eigenvectors.at<float>(0,0));
		float scale = sqrt(mean_curve_pca.eigenvalues.at<float>(0,0));
		Point2f cp = Point2f(mean_curve_pca.mean.at<float>(0),mean_curve_pca.mean.at<float>(1));
        
		cout << "cp " << cp << " angle " << angle << " scale " << scale << "\n";
		
        Mat_<float> Tpca(3,3);
		getRotationMatrix2D(Point2f(0,0), angle * 180.0 / CV_PI, 1.0).copyTo(Tpca.rowRange(0,2));
		Tpca *= (Mat_<float>(3,3) << 1,0,-cp.x, 0,1,-cp.y, 0,0,1);
		//	cout << "rot mat from angle\n"<<getRotationMatrix2D(cp, angle * 180.0 / CV_PI, 1.0) << "\n";
		//	Mat(mean_curve_pca.eigenvectors.t()).copyTo(Tpca.colRange(0, 2));
		
        cout << "trans\n"<<Tpca<<"\n";
		
        transform(mean_curve,mean_curve,Tpca.rowRange(0,2));
		
		Mat show(500,500,CV_8UC3,Scalar::all(0));
		//	drawOpenCurve(show, mean_curve, Scalar(255), 2);
		
//		aligned_curves_by_row.create(aligned_curves.size(),curve_length*2, CV_32FC1);
        assert(aligned_curves_by_row.rows == aligned_curves.size());
        aligned_curves_by_row.setTo(0);
        
		for (int i=0; i<aligned_curves.size(); i++) {
            if (i==first_curve) {
                continue;
            }
//            show.setTo(0);
//			drawOpenCurve(show, mean_curve_orig, Scalar(255,255,255), 1, cv::Point(150,150));
//			drawOpenCurve(show, aligned_curves[i], Scalar(255,0,255), 1, cv::Point(150,150));
            
			vector<Point2f> transformed(mean_curve.size());
            Mat Trans = Find2DRigidTransform(aligned_curves[i], mean_curve_orig);
			transform(aligned_curves[i],transformed,Trans);
            
            aligned_curves[i].clear(); aligned_curves[i].insert(aligned_curves[i].begin(),transformed.begin(),transformed.end());
            
//			drawOpenCurve(show, aligned_curves[i], Scalar(0,0,255), 1, cv::Point(150,150));
//            imshow("aligned",show);
//            waitKey();
            
            Mat transformed_mat = Mat(transformed).t();
            Mat transformed_mat_one_ch = transformed_mat.reshape(1, 1);
            if(!isnan(transformed_mat_one_ch.at<float>(0)))
                transformed_mat_one_ch.copyTo(aligned_curves_by_row.row(i));
        }
		
		Mat curves_mean;
		curves_mean_v.clear(); curves_mean_v.resize(curve_length);
		reduce(aligned_curves_by_row, curves_mean, 0, CV_REDUCE_AVG);
        cout << "curves_mean " << curves_mean << endl;
        Mat curves_mean_t = curves_mean.t();
        Mat curves_mean_t_2ch = curves_mean_t.reshape(2,curve_length);
		curves_mean_t_2ch.copyTo(Mat(curves_mean_v));
        
        original_curves.clear(); original_curves.insert(original_curves.begin(),aligned_curves.begin(),aligned_curves.end());
        buildHPDM();
    }
    
    void saveToFile(const char* filename) {
        FileStorage fs(filename,FileStorage::WRITE);
//        vector<vector<Point_<T> > > original_curves;
//        vector<vector<Point_<T> > > aligned_curves;
//        Mat aligned_curves_by_row;
//        int curve_length;
//        int shape_space_dim;
//        vector<Point2f> curves_mean_v;
//        
//        PCA curves_pca;
//        vector<PCA> subregions_pca;
//        Mat subregion_means;
//        vector<T> subregions_lambda;
//        fs << "original_curves" << original_curves;
//        fs << "aligned_curves" << aligned_curves;

        fs << "curve_length" << curve_length;
        fs << "shape_space_dim" << shape_space_dim;
        fs << "num_pc_for_subregions" << num_pc_for_subregions;
        fs << "aligned_curves_by_row" << aligned_curves_by_row;
        fs << "curves_mean_v" << curves_mean_v;
        fs << "curves_pca_eval" << curves_pca.eigenvalues;
        fs << "curves_pca_evec" << curves_pca.eigenvectors;
        fs << "curves_pca_mean" << curves_pca.mean;
        fs << "subregions_pca" << "[";
        for (int i=0; i<subregions_pca.size(); i++) {
            fs << "{:";
            fs << "eigenvalues" << subregions_pca[i].eigenvalues;
            fs << "eigenvectors" << subregions_pca[i].eigenvectors;
            fs << "mean" << subregions_pca[i].mean;
            fs << "}";
        }
        fs << "]";
        fs << "subregion_means" << subregion_means;
        fs << "subregions_lambda" << subregions_lambda;
        fs.release();
    }
    
    void initializeFromFile(const char* filename) {
        clearDataStruct();
        
        FileStorage fs(filename,FileStorage::READ);
        if(!fs.isOpened()) {
            cerr << "cannot read " << filename << endl;
            return;
        }
//        fs["original_curves"] >> original_curves;
//        fs["aligned_curves" ] >> aligned_curves;
        fs["curve_length" ] >> curve_length;
        fs["shape_space_dim" ] >> shape_space_dim;
        fs["num_pc_for_subregions" ] >> num_pc_for_subregions;
        fs["aligned_curves_by_row"] >> aligned_curves_by_row;
        for (int i=0; i<aligned_curves_by_row.rows; i++) {
            vector<Point_<T> > aligned_c(aligned_curves_by_row.cols / 2);
            Mat(aligned_curves_by_row.row(i).t()).reshape(2,curve_length).copyTo(Mat(aligned_c));
            aligned_curves.push_back(aligned_c);
        }
        fs["curves_mean_v"] >> curves_mean_v;
        fs["curves_pca_eval" ] >> curves_pca.eigenvalues;
        fs["curves_pca_evec" ] >> curves_pca.eigenvectors;
        fs["curves_pca_mean" ] >> curves_pca.mean;
        FileNode features = fs["subregions_pca"];
        FileNodeIterator it = features.begin(), it_end = features.end();
        int idx = 0;
        // iterate through a sequence using FileNodeIterator
        for( ; it != it_end; ++it, idx++ )
        {
            PCA pca;
            (*it)["eigenvalues"] >> pca.eigenvalues;
            (*it)["eigenvectors"] >> pca.eigenvectors;
            (*it)["mean"] >> pca.mean;
            subregions_pca.push_back(pca);
        }
        
        fs["subregion_means" ] >> subregion_means;
        fs["subregions_lambda"] >>  subregions_lambda;
        fs.release();
        
        this->initialized = true;
    }
    
	void initialize(const vector<vector<Point_<T> > >& curves) {
		assert(curves.size() > 0);
		
		original_curves.clear();
		original_curves.insert(original_curves.begin(),curves.begin(),curves.end());
		curve_length = curves[0].size();
		
		align();
		buildHPDM();
		
		this->initialized = true;
	}
	
	int nearestPatch(const vector<Point_<T> >& curve, Mat_<double>& dists) {
		assert(curve.size() == curve_length);
		
		//align
		Mat Trans = Find2DRigidTransform(curve, curves_mean_v);
		vector<Point2f> transformed(curve_length);
		transform(curve,transformed,Trans);
		
		
		//project to shape space
		Mat curve_in_shape_space;
		Mat curve_m = Mat(Mat(transformed).t()).reshape(1,1);
		curves_pca.project(curve_m,curve_in_shape_space);
		
        return nearestPatchWithShapeSpaceCurve(curve_in_shape_space, dists);
    }
    
    int nearestPatchWithShapeSpaceCurve(const Mat& curve_in_shape_space, Mat_<double>& dists) {
		dists.create(1, subregion_means.rows);
		for (int i=0; i<subregion_means.rows; i++) {
			double d = norm(subregion_means.row(i)-curve_in_shape_space);
			double lambda_sq = subregions_lambda[i]; //*subregions_lambda[i];
			double denom = lambda_sq*sqrt(2.0*CV_PI);
			double dist = exp(-d/lambda_sq)/denom;

			dists(i) = dist;
		}
        cv::Point maxpatch; double minVal;
		minMaxLoc(dists, &minVal, 0, 0, &maxpatch);
		return maxpatch.x;
	}
	
    vector<vector<Point_<T> > > getAlignedCurves() { return aligned_curves; };
	int getShapeSpaceDim() const { return shape_space_dim; }
    int getCurveLength() const { return curve_length;}
	int getNumPatches() const { return subregions_pca.size(); }
	Mat getPatchMean(int ptch) const { assert(ptch>=0 && ptch<subregion_means.rows); return subregion_means.row(ptch); }
	vector<Point_<T> > rebuildCurve(const Mat& shape_space_deformation) const { 
		assert(shape_space_deformation.rows == 1 && shape_space_deformation.cols == shape_space_dim);
		
		Mat rebuilt; curves_pca.backProject(shape_space_deformation,rebuilt);
		vector<Point2f> rebuilt_curve(curve_length);
		Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));	
		return rebuilt_curve;
	}
	Mat getPatchPrincipalAxis(int ptch, float& sigma) const {
		assert(ptch>=0 && ptch<subregion_means.rows);
//		Mat pa(1,shape_space_dim,CV_32FC1);
//		Mat().copyTo(pa);
		sigma = sqrt(subregions_pca[ptch].eigenvalues.template at<float>(0));
		return subregions_pca[ptch].eigenvectors.row(0);
	}
    Mat projectToShapeSpace(const vector<Point_<T> >& curve) const {
        assert(curve.size()==curve_length);
        return curves_pca.project(Mat(Mat(curve).t()).reshape(1,1));
    }
    Mat projectToPatchSpace(Mat& in_shape_space, int patch_num) const {
        return subregions_pca[patch_num].project(in_shape_space);
    }
    Mat backprojectFromPatchSpace(Mat& in_patch_space, int patch_num) const {
        return subregions_pca[patch_num].backProject(in_patch_space);
    }
    
    void visualizePatches() const {
        int squaresize = 200;
        int sq_k = cvCeil(sqrt(subregion_means.rows));
		Mat show(squaresize*sq_k,squaresize*sq_k,CV_8UC3,Scalar::all(0));
        
        for (int i=0; i<subregion_means.rows; i++) {
            Mat showroi = show(Range(squaresize*(i/sq_k),squaresize*(i/sq_k + 1)),Range(squaresize*(i%sq_k),squaresize*(i%sq_k + 1)));
            
            Mat rebuilt;
            curves_pca.backProject(subregion_means.row(i),rebuilt);
            vector<Point2f> rebuilt_curve(curve_length);
            Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
            Mat(rebuilt_curve) += Scalar(squaresize/2,squaresize/2);
            drawOpenCurve(showroi, rebuilt_curve, Scalar(0,0,255), 1);
            
            curves_pca.backProject(subregion_means.row(i) + subregions_pca[i].eigenvectors.row(0) * sqrt(subregions_pca[i].eigenvalues.template at<float>(0)),rebuilt);
            Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
            Mat(rebuilt_curve) += Scalar(squaresize/2,squaresize/2);
            drawOpenCurve(showroi, rebuilt_curve, Scalar(0,255,255), 1);
            
            curves_pca.backProject(subregion_means.row(i) - subregions_pca[i].eigenvectors.row(0) * sqrt(subregions_pca[i].eigenvalues.template at<float>(0)),rebuilt);
            Mat(rebuilt).reshape(2,curve_length).copyTo(Mat(rebuilt_curve));
            Mat(rebuilt_curve) += Scalar(squaresize/2,squaresize/2);
            drawOpenCurve(showroi, rebuilt_curve, Scalar(0,255,255), 1);
            
            stringstream ss; ss << "patch " << i;
            putText(showroi, ss.str(), cv::Point(10,10), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(255));
            
        }
        imshow("subregions means",show);
        waitKey(0);
    }

};