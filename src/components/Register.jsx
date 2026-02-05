import React, { useState } from 'react';

const Register = ({ onRegister, onSwitchToLogin }) => {
    const [formData, setFormData] = useState({
        firstName: '',
        lastName: '',
        email: '',
        password: '',
        confirmPassword: '',
        organization: ''
    });

    const handleChange = (e) => {
        setFormData({
            ...formData,
            [e.target.name]: e.target.value
        });
    };

    const handleSubmit = (e) => {
        e.preventDefault();
        if (formData.password !== formData.confirmPassword) {
            alert("Passwords do not match!");
            return;
        }
        console.log('Register attempt:', formData);
        if (onRegister) onRegister(formData);
    };

    return (
        <div style={styles.container}>
            <div style={styles.card}>
                <div style={styles.header}>
                    <h2 style={styles.title}>Create Account</h2>
                    <p style={styles.subtitle}>Join MoleculeAI to start discovering</p>
                </div>

                <form onSubmit={handleSubmit} style={styles.form}>
                    <div style={styles.row}>
                        <div style={styles.inputGroup}>
                            <label style={styles.label}>First Name</label>
                            <input
                                type="text"
                                name="firstName"
                                value={formData.firstName}
                                onChange={handleChange}
                                placeholder="First Name"
                                style={styles.input}
                                required
                            />
                        </div>
                        <div style={styles.inputGroup}>
                            <label style={styles.label}>Last Name</label>
                            <input
                                type="text"
                                name="lastName"
                                value={formData.lastName}
                                onChange={handleChange}
                                placeholder="Last Name"
                                style={styles.input}
                                required
                            />
                        </div>
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Email Address</label>
                        <input
                            type="email"
                            name="email"
                            value={formData.email}
                            onChange={handleChange}
                            placeholder="Enter your email"
                            style={styles.input}
                            required
                        />
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Password</label>
                        <input
                            type="password"
                            name="password"
                            value={formData.password}
                            onChange={handleChange}
                            placeholder="Create a password"
                            style={styles.input}
                            required
                        />
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Confirm Password</label>
                        <input
                            type="password"
                            name="confirmPassword"
                            value={formData.confirmPassword}
                            onChange={handleChange}
                            placeholder="Confirm your password"
                            style={styles.input}
                            required
                        />
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>
                            Organization <span style={styles.optional}>(Optional)</span>
                        </label>
                        <input
                            type="text"
                            name="organization"
                            value={formData.organization}
                            onChange={handleChange}
                            placeholder="Company or University"
                            style={styles.input}
                        />
                    </div>

                    <button type="submit" style={styles.button}>
                        Create Account
                    </button>
                </form>

                <div style={styles.switchSection}>
                    <span style={styles.switchText}>Already have an account?</span>
                    <button
                        type="button"
                        style={styles.switchButton}
                        onClick={onSwitchToLogin}
                    >
                        Sign In
                    </button>
                </div>
            </div>
        </div>
    );
};

const styles = {
    container: {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        minHeight: '700px',
        padding: '20px',
        animation: 'fadeIn 0.5s ease-out',
    },
    card: {
        width: '100%',
        maxWidth: '500px',
        backgroundColor: '#ffffff',
        borderRadius: '16px',
        boxShadow: '0 4px 20px rgba(0, 0, 0, 0.05)',
        padding: '40px',
        border: '1px solid #f3f4f6',
    },
    header: {
        textAlign: 'center',
        marginBottom: '32px',
    },
    title: {
        fontSize: '28px',
        fontWeight: '700',
        color: '#111827',
        margin: '0 0 8px 0',
        fontFamily: "'DM Sans', sans-serif",
    },
    subtitle: {
        fontSize: '15px',
        color: '#6b7280',
        margin: 0,
        lineHeight: '1.5',
        fontFamily: "'DM Sans', sans-serif",
    },
    form: {
        display: 'flex',
        flexDirection: 'column',
        gap: '20px',
    },
    row: {
        display: 'flex',
        gap: '16px',
    },
    inputGroup: {
        display: 'flex',
        flexDirection: 'column',
        gap: '8px',
        flex: 1,
    },
    label: {
        fontSize: '14px',
        fontWeight: '600',
        color: '#374151',
        fontFamily: "'DM Sans', sans-serif",
    },
    optional: {
        fontWeight: '400',
        color: '#9ca3af',
        fontSize: '13px',
    },
    input: {
        padding: '12px 16px',
        borderRadius: '8px',
        border: '1px solid #e5e7eb',
        fontSize: '15px',
        backgroundColor: '#f9fafb',
        transition: 'all 0.2s ease',
        outline: 'none',
        fontFamily: "'DM Sans', sans-serif', width: '100%'",
    },
    button: {
        marginTop: '12px',
        padding: '14px',
        backgroundColor: '#16a34a',
        color: '#ffffff',
        border: 'none',
        borderRadius: '8px',
        fontSize: '16px',
        fontWeight: '600',
        cursor: 'pointer',
        transition: 'background-color 0.2s ease',
        fontFamily: "'DM Sans', sans-serif",
    },
    switchSection: {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        gap: '8px',
        marginTop: '24px',
        paddingTop: '20px',
        borderTop: '1px solid #f3f4f6',
    },
    switchText: {
        fontSize: '14px',
        color: '#6b7280',
        fontFamily: "'DM Sans', sans-serif",
    },
    switchButton: {
        fontSize: '14px',
        fontWeight: '600',
        color: '#16a34a',
        background: 'none',
        border: 'none',
        cursor: 'pointer',
        fontFamily: "'DM Sans', sans-serif",
        padding: 0,
        transition: 'color 0.2s ease',
    }
};

export default Register;
