use ark_bls12_381::Bls12_381;
use ark_ec::pairing::{Pairing, self};
use ark_ff::UniformRand;
use ark_ip_proofs::applications::poly_commit::{
    transparent::UnivariatePolynomialCommitment as TransparentIPA,
    UnivariatePolynomialCommitment as IPA, KZG,
};
use ark_poly::polynomial::{
    univariate::DensePolynomial as UnivariatePolynomial, DenseUVPolynomial, Polynomial,
};

use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use csv::Writer;

use blake2::Blake2b;
use digest::Digest;
use std::{
    io::stdout,
    time::{Duration, Instant}, vec,
    mem::{size_of, size_of_val},
};

use ark_ff::Fp;
use ark_bls12_381::FrConfig;
use ark_ff::MontBackend;
use ark_bls12_381::Config;


fn main() {
    // // A test for Blake2b
    // let mut hasher = Blake2b::new();
    // // write input message
    // hasher.update(b"hello world");
    // // read hash digest and consume hasher
    // let res = hasher.finalize();
    // println!("GenericArray hash: {:?}", res);


    let mut args: Vec<String> = std::env::args().collect();
    if args.last().unwrap() == "--bench" {
        args.pop();
    }
    let (num_trials, num_data_points): (usize, usize) = if args.len() < 2
        || args[1] == "-h"
        || args[1] == "--help"
    {
        println!("Usage: ``cargo bench --bench poly_commit -- <num_trials> <num_data_points/degree_bound>``");
        return;
    } else {
        (
            String::from(args[1].clone())
                .parse()
                .expect("<num_trials> should be integer"),
            String::from(args[2].clone())
                .parse()
                .expect("<num_data_points> should be integer"),
        )
    };

    let mut csv_writer = Writer::from_writer(stdout());
    csv_writer
        .write_record(&["testnum", "scheme", "function", "degree/size (bytes)", "time (ms)"])
        .unwrap();
    csv_writer.flush().unwrap();

    let mut start;
    let mut time;

    println!("field_size: 32 bytes");
    println!("G1_size: 144 bytes");
    println!("G2_size: 288 bytes");
    println!("GT_size: 576 bytes");

    // degree from (4^1 - 1) to (4^num_data_points - 1)
    // 4 guarantees that sqrt is power-of-two
    for degree in (13..num_data_points).map(|i| 4_usize.pow((i + 1) as u32) - 1) {
        // Benchmark nameKZG
        {

        let degree_sqrt_f64 = (((degree + 1) * 64) as f64).sqrt();
        let degree_sqrt = (degree_sqrt_f64 as usize) - 1;
        println!("degree_sqrt: {}", degree_sqrt);
        let poly_num = (degree + 1) / (degree_sqrt + 1);
        println!("poly_num: {}", poly_num);

        let mut rng = StdRng::seed_from_u64(0u64);

        start = Instant::now();
        let (g_alpha_powers, v_srs) = KZG::<Bls12_381>::setup(&mut rng, degree_sqrt).unwrap();
        time = start.elapsed().as_millis();

        // record element size
        // let v_srs_size = size_of_val(&v_srs);
        // println!("v_srs_size: {}", v_srs_size);
        let field_size = size_of_val(&<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        // println!("field_size: {}", field_size);
        let G1_size = size_of_val(&v_srs.g);
        // println!("G1_size: {}", G1_size);
        let G2_size = size_of_val(&v_srs.h);
        // println!("G2_size: {}", G2_size);

        let srs_size = (g_alpha_powers.len() + 2) * G1_size + 2 * G2_size;

        csv_writer
            .write_record(&[
                1.to_string(),
                "namekzg".to_string(),
                "setup".to_string(),
                srs_size.to_string(),
                time.to_string(),
            ])
            .unwrap();
        csv_writer.flush().unwrap();

        for i in 1..num_trials + 1 {

            let mut polynomial_vec: Vec<ark_poly::univariate::DensePolynomial<Fp<MontBackend<FrConfig, 4>, 4>>> = Vec::with_capacity(poly_num);

            for j in 0..(poly_num) {
                let polynomial = UnivariatePolynomial::rand(degree_sqrt, &mut rng);
                // polynomial_coeff_vec
                polynomial_vec.push(polynomial);
            }

            let point = <Bls12_381 as Pairing>::ScalarField::rand(&mut rng);

            // let mut eval_vec: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>> = Vec::with_capacity(poly_num);
            let mut final_eva = polynomial_vec[0].evaluate(&point);
            for j in 1..(poly_num) {
                final_eva += polynomial_vec[j].evaluate(&point);
            }

            // Commit

            let mut com_vec: Vec<<Bls12_381 as Pairing>::G1> = Vec::with_capacity(poly_num);
            let mut field_vec: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>> = Vec::with_capacity(poly_num);

            start = Instant::now();

            for j in 0..(poly_num) {
                com_vec.push(KZG::<Bls12_381>::commit(&g_alpha_powers, &polynomial_vec[j]).unwrap());
                let mut hasher = Blake2b::new();
                // write input message
                hasher.update(com_vec[j].to_string());
                // read hash digest and consume hasher
                let res = hasher.finalize();
                let res_field: ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4> = (0..16).into_iter().fold(0u128, |sum, next |sum * 256 + res[next] as u128).into();
                field_vec.push(res_field);
            }

            let com_poly = UnivariatePolynomial::from_coefficients_vec(field_vec); 
            let com_poly_com = KZG::<Bls12_381>::commit(&g_alpha_powers, &com_poly).unwrap();

            time = start.elapsed().as_millis();

            csv_writer
                .write_record(&[
                    i.to_string(),
                    "namekzg".to_string(),
                    "commit".to_string(),
                    size_of_val(&com_poly_com).to_string(),
                    time.to_string(),
                ])
                .unwrap();

            // Open

            let mut batch_polynomial = &polynomial_vec[0] + &polynomial_vec[1];
            let mut batch_com: <Bls12_381 as Pairing>::G1 = com_vec[0] + com_vec[1];

            start = Instant::now();

            for j in 2..(poly_num) {
                batch_polynomial += &polynomial_vec[j];
                batch_com += com_vec[j];
            }

            let proof = (KZG::<Bls12_381>::open(&g_alpha_powers, &batch_polynomial, &point).unwrap());

            time = start.elapsed().as_millis();

            let proof_size = poly_num * (G1_size + field_size) + G1_size;

            csv_writer
                .write_record(&[
                    i.to_string(),
                    "namekzg".to_string(),
                    "open".to_string(),
                    proof_size.to_string(),
                    time.to_string(),
                ])
                .unwrap();

            // Verify
            std::thread::sleep(Duration::from_millis(5000));

            start = Instant::now();
            for _ in 0..50 {
                let mut field_vec: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>> = Vec::with_capacity(poly_num);
                for j in 0..(poly_num) {
                    let mut hasher = Blake2b::new();
                    // write input message
                    hasher.update(com_vec[j].to_string());
                    // read hash digest and consume hasher
                    let res = hasher.finalize();
                    let res_field: ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4> = (0..16).into_iter().fold(0u128, |sum, next |sum * 256 + res[next] as u128).into();
                    field_vec.push(res_field);
                }

                let com_poly = UnivariatePolynomial::from_coefficients_vec(field_vec); 
                let com_poly_com = KZG::<Bls12_381>::commit(&g_alpha_powers, &com_poly).unwrap();

                let is_valid =
                    KZG::<Bls12_381>::verify(&v_srs, &batch_com, &point, &final_eva, &proof).unwrap();
                assert!(is_valid);
            }
            time = start.elapsed().as_millis() / 50;
            csv_writer
                .write_record(&[
                    i.to_string(),
                    "namekzg".to_string(),
                    "verify".to_string(),
                    time.to_string(),
                    "ms".to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            }
        }
        
        // Benchmark KZG
        {
            let mut rng = StdRng::seed_from_u64(0u64);
            println!("degree: {}", degree);
            start = Instant::now();
            let (g_alpha_powers, v_srs) = KZG::<Bls12_381>::setup(&mut rng, degree).unwrap();
            time = start.elapsed().as_millis();

            // record element size
            // let v_srs_size = size_of_val(&v_srs);
            // println!("v_srs_size: {}", v_srs_size);
            let field_size = size_of_val(&<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            // println!("field_size: {}", field_size);
            let G1_size = size_of_val(&v_srs.g);
            // println!("G1_size: {}", G1_size);
            let G2_size = size_of_val(&v_srs.h);
            // println!("G2_size: {}", G2_size);

            let srs_size = (g_alpha_powers.len() + 2) * G1_size + 2 * G2_size;
            // println!("srs_size: {}", srs_size);
            csv_writer
                .write_record(&[
                    1.to_string(),
                    "kzg".to_string(),
                    "setup".to_string(),
                    srs_size.to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            for i in 1..num_trials + 1 {
                let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
                let point = <Bls12_381 as Pairing>::ScalarField::rand(&mut rng);
                let eval = polynomial.evaluate(&point);

                // Commit
                start = Instant::now();
                let com = KZG::<Bls12_381>::commit(&g_alpha_powers, &polynomial).unwrap();

                time = start.elapsed().as_millis();
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "commit".to_string(),
                        size_of_val(&com).to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Open
                start = Instant::now();
                let proof = KZG::<Bls12_381>::open(&g_alpha_powers, &polynomial, &point).unwrap();
                time = start.elapsed().as_millis();
                let proof_size = size_of_val(&proof);
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "open".to_string(),
                        proof_size.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Verify
                std::thread::sleep(Duration::from_millis(5000));
                start = Instant::now();
                for _ in 0..50 {
                    let is_valid =
                        KZG::<Bls12_381>::verify(&v_srs, &com, &point, &eval, &proof).unwrap();
                    assert!(is_valid);
                }
                time = start.elapsed().as_millis() / 50;
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "kzg".to_string(),
                        "verify".to_string(),
                        time.to_string(),
                        "ms".to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
            }
        }

        // Benchmark IPA
        {
            let mut rng = StdRng::seed_from_u64(0u64);
            println!("IPA_degree: {}", degree);
            start = Instant::now();
            let srs = IPA::<Bls12_381, Blake2b>::setup(&mut rng, degree).unwrap();
            let v_srs = srs.0.get_verifier_key();
            time = start.elapsed().as_millis();

            let field_size = size_of_val(&<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            // println!("field_size: {}", field_size);
            // println!("kzg_srs_size: {}", srs.1.len());
            let com_key = srs.0.get_commitment_keys();
            // println!("com_key_size: {}", com_key.0.len());
            let G1_size = size_of_val(&srs.1[0]);
            let srs_G1_size = size_of_val(&srs.1[0]) * (srs.1.len() + 2);
            // println!("G1_size: {}", G1_size);
            let G2_size = size_of_val(&com_key.0[0]);
            let srs_G2_size = size_of_val(&com_key.0[0]) * (com_key.0.len() + 2);
            // println!("G2_size: {}", G2_size);

            csv_writer
                .write_record(&[
                    1.to_string(),
                    "ipa".to_string(),
                    "setup".to_string(),
                    (srs_G1_size + srs_G2_size).to_string(),
                    time.to_string(),
                ])
                .unwrap();
            csv_writer.flush().unwrap();
            for i in 1..num_trials + 1 {
                let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
                let point = <Bls12_381 as Pairing>::ScalarField::rand(&mut rng);
                let eval = polynomial.evaluate(&point);

                // Commit
                start = Instant::now();
                let (com, prover_aux) =
                    IPA::<Bls12_381, Blake2b>::commit(&srs, &polynomial).unwrap();
                time = start.elapsed().as_millis();

                let GT_size = size_of_val(&com);
                // println!("GT_size: {}", GT_size);

                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "commit".to_string(),
                        GT_size.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Open
                start = Instant::now();
                let proof = IPA::<Bls12_381, Blake2b>::open(&srs, &polynomial, &prover_aux, &point)
                    .unwrap();
                time = start.elapsed().as_millis();

                // proof_size = 2 logm GT + 1 G1 + 2 G2
                // kzg_size = 1 G1
                
                let proof_size = GT_size * 2 * (((com_key.0.len()/2) as f64).log2() as usize) + 2 * G1_size + 2 * G2_size;

                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "open".to_string(),
                        proof_size.to_string(),
                        time.to_string(),
                    ])
                    .unwrap();

                // Verify
                std::thread::sleep(Duration::from_millis(5000));
                start = Instant::now();
                for _ in 0..50 {
                    let is_valid = IPA::<Bls12_381, Blake2b>::verify(
                        &v_srs, degree, &com, &point, &eval, &proof,
                    )
                    .unwrap();
                    assert!(is_valid);
                }
                time = start.elapsed().as_millis() / 50;
                csv_writer
                    .write_record(&[
                        i.to_string(),
                        "ipa".to_string(),
                        "verify".to_string(),
                        time.to_string(),
                        "ms".to_string(),
                    ])
                    .unwrap();
                csv_writer.flush().unwrap();
            }
        }

            // // Benchmark transparent IPA
            // {
            //     let mut rng = StdRng::seed_from_u64(0u64);
            //     start = Instant::now();
            //     let ck = TransparentIPA::<Bls12_381, Blake2b>::setup(&mut rng, degree).unwrap();
            //     time = start.elapsed().as_millis();
            //     csv_writer
            //         .write_record(&[
            //             1.to_string(),
            //             "transparent_ipa".to_string(),
            //             "setup".to_string(),
            //             degree.to_string(),
            //             time.to_string(),
            //         ])
            //         .unwrap();
            //     csv_writer.flush().unwrap();
            //     for i in 1..num_trials + 1 {
            //         let polynomial = UnivariatePolynomial::rand(degree, &mut rng);
            //         let point = <Bls12_381 as Pairing>::ScalarField::rand(&mut rng);
            //         let eval = polynomial.evaluate(&point);

            //         // Commit
            //         start = Instant::now();
            //         let (com, prover_aux) =
            //             TransparentIPA::<Bls12_381, Blake2b>::commit(&ck, &polynomial).unwrap();
            //         time = start.elapsed().as_millis();
            //         csv_writer
            //             .write_record(&[
            //                 i.to_string(),
            //                 "transparent_ipa".to_string(),
            //                 "commit".to_string(),
            //                 degree.to_string(),
            //                 time.to_string(),
            //             ])
            //             .unwrap();

            //         // Open
            //         start = Instant::now();
            //         let proof = TransparentIPA::<Bls12_381, Blake2b>::open(
            //             &ck,
            //             &polynomial,
            //             &prover_aux,
            //             &point,
            //         )
            //         .unwrap();
            //         time = start.elapsed().as_millis();
            //         csv_writer
            //             .write_record(&[
            //                 i.to_string(),
            //                 "transparent_ipa".to_string(),
            //                 "open".to_string(),
            //                 degree.to_string(),
            //                 time.to_string(),
            //             ])
            //             .unwrap();

            //         // Verify
            //         std::thread::sleep(Duration::from_millis(5000));
            //         start = Instant::now();
            //         for _ in 0..50 {
            //             let is_valid = TransparentIPA::<Bls12_381, Blake2b>::verify(
            //                 &ck, &com, &point, &eval, &proof,
            //             )
            //             .unwrap();
            //             assert!(is_valid);
            //         }
            //         time = start.elapsed().as_millis() / 50;
            //         csv_writer
            //             .write_record(&[
            //                 i.to_string(),
            //                 "transparent_ipa".to_string(),
            //                 "verify".to_string(),
            //                 degree.to_string(),
            //                 time.to_string(),
            //             ])
            //             .unwrap();
            //         csv_writer.flush().unwrap();
            //     }
            // }
        // }
    }
}
